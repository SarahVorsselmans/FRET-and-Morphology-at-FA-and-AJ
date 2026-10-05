function [sigPairs, ranovaTbl, rmModel, wideTbl] = runRepeatedMeasuresANOVA(tbl, metricVec, options)
%RUNREPEATEDMEASURESANOVA
% Repeated-measures ANOVA on bioRep means (fully paired design required)
%
% REQUIREMENT:
%   Each bioRep must have ALL conditions (you said you'll pre-filter)
%
% -------------------------------------------------------------------------

arguments
    tbl table
    metricVec (:,1) double

    options.CombinedVec           = {}
    options.BioRepCol string      = "expCat"
    options.ConditionLevels cell  = {}
    options.Correction string     = "fdr"
    options.Alpha (1,1) double    = 0.05
    options.Verbose (1,1) logical = true
end

alphaVal = options.Alpha;

%% --- Extract vectors ---
bioVec  = string(tbl.(options.BioRepCol));
condVec = string(options.CombinedVec);

if isempty(options.ConditionLevels)
    condLevels = cellstr(unique(condVec(~isnan(metricVec)), 'stable'));
else
    condLevels = options.ConditionLevels;
end

nCond = numel(condLevels);

%% --- Build long table ---
validMask = ~isnan(metricVec) & ismember(condVec, condLevels);

longTbl = table();
longTbl.metric    = metricVec(validMask);
longTbl.condition = categorical(condVec(validMask), condLevels);
longTbl.bioRep    = categorical(bioVec(validMask));

%% --- Collapse to bioRep means (required for rmANOVA) ---
bioLevels = categories(longTbl.bioRep);

rows = {};

validBio = true(numel(bioLevels),1);

% --- Step 1: check which bioReps are fully paired ---
for b = 1:numel(bioLevels)
    for k = 1:nCond
        mask = longTbl.bioRep == bioLevels{b} & ...
               longTbl.condition == condLevels{k};

        vals = longTbl.metric(mask);

        if isempty(vals)
            validBio(b) = false;
            break;
        end
    end
end

bioLevels = bioLevels(validBio);

if isempty(bioLevels)
    error('runRepeatedMeasuresANOVA: no fully paired bioReps found.');
end

% --- Step 2: build rows only from valid bioReps ---
for b = 1:numel(bioLevels)
    for k = 1:nCond
        mask = longTbl.bioRep == bioLevels{b} & ...
               longTbl.condition == condLevels{k};

        vals = longTbl.metric(mask);

        rows(end+1,:) = { ...
            mean(vals,'omitnan'), ...
            condLevels{k}, ...
            bioLevels{b} ...
        }; %#ok<AGROW>
    end
end

% --- Step 3: convert safely ---
tmp = rows;

aggTbl = table();
aggTbl.metric    = cell2mat(tmp(:,1));
aggTbl.condition = categorical(tmp(:,2), condLevels);
aggTbl.bioRep    = categorical(tmp(:,3));

%% --- Convert to wide format ---
wideTbl = unstack(aggTbl, 'metric', 'condition');

% Ensure ordering
wideTbl = sortrows(wideTbl, 'bioRep');

%% --- Sanitize column names for fitrm compatibility ---
% fitrm requires valid MATLAB identifiers as column names (no spaces, special chars)
% Map original condLevels -> safe names, then rename wideTbl columns accordingly

safeNames = matlab.lang.makeValidName(condLevels);
safeNames = matlab.lang.makeUniqueStrings(safeNames);  % guard against collisions

% Build a rename map from whatever unstack produced -> safeNames
% unstack uses the category labels as column names (may already be mangled by MATLAB)
% Safer: just force-assign the variable names directly (skip bioRep col)
bioRepColIdx = strcmp(wideTbl.Properties.VariableNames, 'bioRep');
dataColIdx   = ~bioRepColIdx;
wideTbl.Properties.VariableNames(dataColIdx) = safeNames;

%% --- Define repeated measures model ---
% Use safeNames for the formula and WithinDesign; keep condLevels for output
withinDesign = table(categorical(condLevels'), 'VariableNames', {'Condition'});

formula = sprintf('%s-%s ~ 1', safeNames{1}, safeNames{end});

rmModel = fitrm(wideTbl, formula, 'WithinDesign', withinDesign);

% %% --- Define repeated measures model ---
% measNames = condLevels;
% 
% withinDesign = table(categorical(measNames'), 'VariableNames', {'Condition'});
% 
% formula = sprintf('%s-%s ~ 1', measNames{1}, measNames{end});
% 
% rmModel = fitrm(wideTbl, formula, 'WithinDesign', withinDesign);

%% --- rmANOVA ---
ranovaTbl = ranova(rmModel);

try
    sphTbl = mauchly(rmModel);
    sphericityP = sphTbl.pValue(end);
    if options.Verbose
        fprintf('  Mauchly: W=%.4f, p=%.4g\n', sphTbl.W(end), sphericityP);
    end
    if sphericityP < 0.05
        epsTbl = epsilon(rmModel);
        gg  = epsTbl.GreenhouseGeisser(end);
        effectRow = contains(ranovaTbl.Properties.RowNames, ':') & ...
                    ~contains(ranovaTbl.Properties.RowNames, 'Intercept');
        errorRow  = ~contains(ranovaTbl.Properties.RowNames, ':') & ...
                     contains(ranovaTbl.Properties.RowNames, 'Error');
        F   = ranovaTbl.F(effectRow);
        df1 = ranovaTbl.DF(effectRow) * gg;
        df2 = ranovaTbl.DF(errorRow)  * gg;
        pGG = 1 - fcdf(F, df1, df2);
        ranovaTbl.pValueGG = NaN(height(ranovaTbl), 1);
        ranovaTbl.pValueGG(effectRow) = pGG;
        if options.Verbose
            fprintf('  Sphericity violated — GG (eps=%.4f): F=%.4f, p=%.4g\n', gg, F, pGG);
        end
    end
catch ME
    warning('%s', sprintf('Sphericity test failed: %s', ME.message));
end

if options.Verbose
    fprintf('\n=== Repeated Measures ANOVA ===\n');
    disp(ranovaTbl);
end

%% --- Pairwise comparisons ---
mc = multcompare(rmModel, 'Condition', 'ComparisonType', options.Correction);
% Then use mc.pValue directly — no further multiplication
adjP = mc.pValue;

comboA = cellstr(string(mc.Condition_1));
comboB = cellstr(string(mc.Condition_2));
rawP   = adjP;

nPairs = nchoosek(nCond, 2);

% %% --- Multiple comparison correction ---
% switch lower(options.Correction)
%     case 'fdr'
%         adjP = mafdr(rawP, 'BHFDR', true);
%     case 'bonferroni'
%         adjP = min(rawP * nPairs, 1);
%     otherwise
%         adjP = rawP;
% end

%% --- Build outputs ---
% sigPairs = {};
% sigLabel = cell(nPairs,1);
sigLabel = arrayfun(@(p) pToStars(p, alphaVal), adjP, 'UniformOutput', false);
% 
% for i = 1:nPairs
%     sigLabel{i} = pToStars(adjP(i), alphaVal);
%     if adjP(i) < alphaVal
%         sigPairs{end+1} = {comboA{i}, comboB{i}, sigLabel{i}}; %#ok<AGROW>
%     end
% end
%% --- Unique pairs using allPairs ---
pairs = allPairs(condLevels);

sigPairs = {};

for i = 1:numel(pairs)
    A = pairs{i}{1};
    B = pairs{i}{2};

    % Find corresponding comparison (either direction)
    idx = (strcmp(comboA, A) & strcmp(comboB, B)) | ...
          (strcmp(comboA, B) & strcmp(comboB, A));

    if any(idx)
        pval = adjP(idx);
        pval = pval(1);  % same value both directions

        label = pToStars(pval, alphaVal);

        if pval < alphaVal
            sigPairs{end+1} = {A, B, label}; %#ok<AGROW>
        end
    end
end

%% --- Verbose output ---
if options.Verbose
    fprintf('\n--- Pairwise comparisons (%s correction) ---\n', options.Correction);
    fprintf('%-20s %-20s %8s %8s %6s\n','A','B','pRaw','pAdj','Sig');
    fprintf('%s\n', repmat('-',1,70));
    for i = 1:nPairs
        fprintf('%-20s %-20s %8.4g %8.4g %s\n', ...
            comboA{i}, comboB{i}, rawP(i), adjP(i), sigLabel{i});
    end
end

end

% =========================================================================
% LOCAL HELPERS
% =========================================================================
function pairs = allPairs(keys)
    n = numel(keys);
    pairs = {};
    for i = 1:n-1
        for j = i+1:n
            pairs{end+1} = {keys{i}, keys{j}}; %#ok<AGROW>
        end
    end
end

function out = pToStars(p, alpha)
%PTOSTARS  Convert p-value to significance stars relative to alpha.
    if nargin < 2, alpha = 0.05; end
    if isnan(p),                          out = 'n/a';
    elseif p < alpha * 0.001 / 0.05,     out = '***';
    elseif p < alpha * 0.01  / 0.05,     out = '**';
    elseif p < alpha,                     out = '*';
    else,                                 out = 'ns';
    end
end