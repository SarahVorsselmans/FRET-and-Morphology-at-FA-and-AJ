function [sigPairs, statsTbl, workTbl] = runContrastTTests(tbl, metricVec, options)
%RUNCONTRASTTESTS  Paired or independent t-tests on specified contrasts
%
% Tests only the user-specified condition pairs (contrasts), with optional
% biological replicate averaging and multiple comparison correction.
%
% -------------------------------------------------------------------------
% OPTIONS (Name-Value)
%   'CombinedVec'      : condition vector (e.g. 'CT_rest')
%   'BioRepCol'        : column for biological replicate ID
%   'Contrasts'        : cell array of pairs, e.g.
%                          {{'CT_rest','CT_active'}, {'HF_rest','HF_active'}}
%   'MeanBioRep'       : true  → collapse to bioRep means before testing
%                        false → use all observations (no bioRep structure)
%   'TestType'         : 'independent' | 'paired'
%                        'paired' requires MeanBioRep = true (one value per
%                        bioRep per condition)
%   'Correction'       : 'fdr' | 'bonferroni' | 'none'
%   'Alpha'            : significance threshold (default 0.05)
%   'Verbose'          : print output
% -------------------------------------------------------------------------

arguments
    tbl       table
    metricVec (:,1) double

    options.CombinedVec                    = {}
    options.BioRepCol     string           = "expCat"
    options.Contrasts     cell             = {}
    options.MeanBioRep    (1,1) logical    = true
    options.TestType      string           = "independent"
    options.Correction    string           = "fdr"
    options.Alpha         (1,1) double     = 0.05
    options.Verbose       (1,1) logical    = true
end

alphaVal  = options.Alpha;
testType  = lower(options.TestType);
contrasts = options.Contrasts;

if isempty(contrasts)
    error('runContrastTTests: ''Contrasts'' must be provided.');
end

%% --- Validate TestType ---
if ~ismember(testType, {'independent','paired'})
    error('runContrastTTests: TestType must be ''independent'' or ''paired''.');
end
if strcmp(testType, 'paired') && ~options.MeanBioRep
    error('runContrastTTests: paired t-test requires MeanBioRep = true.');
end

%% --- Collect all condition labels involved ---
allConds = unique([contrasts{:}]);   % flat list of unique condition names

%% --- Build condition + bioRep vectors ---
condVec = string(options.CombinedVec);
bioVec  = string(tbl.(options.BioRepCol));

validMask = ~isnan(metricVec) & ismember(condVec, allConds);

workTbl           = table();
workTbl.metric    = metricVec(validMask);
workTbl.condition = categorical(condVec(validMask), allConds);
workTbl.bioRep    = categorical(bioVec(validMask));

%% --- Optional: collapse to bioRep means ---
if options.MeanBioRep
    bioLevels       = categories(workTbl.bioRep);
    condLevels      = categories(workTbl.condition);

    rows = {};
    for b = 1:numel(bioLevels)
        for k = 1:numel(condLevels)
            mask = workTbl.bioRep  == bioLevels{b} & ...
                   workTbl.condition == condLevels{k};
            vals = workTbl.metric(mask);
            if isempty(vals), continue; end
            rows{end+1} = {mean(vals,'omitnan'), condLevels{k}, bioLevels{b}}; %#ok<AGROW>
        end
    end

    if isempty(rows)
        error('runContrastTTests: no data left after bioRep aggregation.');
    end
    tmp = vertcat(rows{:});
    workTbl           = table();
    workTbl.metric    = cell2mat(tmp(:,1));
    workTbl.condition = categorical(tmp(:,2), allConds);
    workTbl.bioRep    = categorical(tmp(:,3));

    if options.Verbose
        fprintf('After bioRep averaging: %d rows\n', height(workTbl));
    end
end

%% --- Run t-tests on each contrast ---
nContrasts = numel(contrasts);
comboA     = cell(nContrasts, 1);
comboB     = cell(nContrasts, 1);
rawP       = nan(nContrasts, 1);
nA_vec     = nan(nContrasts, 1);
nB_vec     = nan(nContrasts, 1);

for i = 1:nContrasts
    condA = contrasts{i}{1};
    condB = contrasts{i}{2};

    comboA{i} = condA;
    comboB{i} = condB;

    valsA = workTbl.metric(workTbl.condition == condA);
    valsB = workTbl.metric(workTbl.condition == condB);

    nA_vec(i) = numel(valsA);
    nB_vec(i) = numel(valsB);

    if numel(valsA) < 2 || numel(valsB) < 2
        warning('runContrastTTests: contrast %s vs %s has too few observations; p set to NaN.', ...
            condA, condB);
        rawP(i) = NaN;
        continue
    end

    switch testType
        case 'independent'
            [~, rawP(i)] = ttest2(valsA, valsB);

        case 'paired'
            % Match by bioRep — only shared bioReps are used
            biorepsA = workTbl.bioRep(workTbl.condition == condA);
            biorepsB = workTbl.bioRep(workTbl.condition == condB);
            sharedBR = intersect(categories(biorepsA), categories(biorepsB));

            if numel(sharedBR) < 2
                warning('runContrastTTests: contrast %s vs %s has fewer than 2 shared bioReps; p set to NaN.', ...
                    condA, condB);
                rawP(i) = NaN;
                continue
            end

            pairedA = arrayfun(@(b) ...
                workTbl.metric(workTbl.condition == condA & workTbl.bioRep == b), ...
                categorical(sharedBR), 'UniformOutput', false);
            pairedB = arrayfun(@(b) ...
                workTbl.metric(workTbl.condition == condB & workTbl.bioRep == b), ...
                categorical(sharedBR), 'UniformOutput', false);

            [~, rawP(i)] = ttest(cell2mat(pairedA), cell2mat(pairedB));
    end
end

%% --- Multiple comparison correction ---
validIdx = ~isnan(rawP);
adjP     = nan(nContrasts, 1);

switch lower(options.Correction)
    case 'fdr'
        adjP(validIdx) = mafdr(rawP(validIdx), 'BHFDR', true);
    case 'bonferroni'
        adjP(validIdx) = min(rawP(validIdx) * sum(validIdx), 1);
    otherwise
        adjP(validIdx) = rawP(validIdx);
end

%% --- Build outputs ---
sigPairs = {};
sigLabel = cell(nContrasts, 1);

for i = 1:nContrasts
    sigLabel{i} = pToStars(adjP(i), alphaVal);
    if ~isnan(adjP(i)) && adjP(i) < alphaVal
        sigPairs{end+1} = {comboA{i}, comboB{i}, sigLabel{i}}; %#ok<AGROW>
    end
end

statsTbl = table(comboA, comboB, nA_vec, nB_vec, rawP, adjP, sigLabel, ...
    'VariableNames', {'ComboA','ComboB','nA','nB','pRaw','pAdj','Sig'});

%% --- Verbose output ---
if options.Verbose
    fprintf('\n=== Contrast t-tests (%s, %s correction) ===\n', ...
        testType, options.Correction);
    fprintf('%-20s %-20s %5s %5s %8s %8s %6s\n', ...
        'A','B','nA','nB','pRaw','pAdj','Sig');
    fprintf('%s\n', repmat('-',1,76));
    for i = 1:nContrasts
        fprintf('%-20s %-20s %5d %5d %8.4g %8.4g %s\n', ...
            comboA{i}, comboB{i}, nA_vec(i), nB_vec(i), ...
            rawP(i), adjP(i), sigLabel{i});
    end
end

end


% =========================================================================
% LOCAL HELPERS
% =========================================================================

function out = pToStars(p, alpha)
    if nargin < 2, alpha = 0.05; end
    if isnan(p),                    out = 'n/a';
    elseif p < alpha * 0.001/0.05,  out = '***';
    elseif p < alpha * 0.01/0.05,   out = '**';
    elseif p < alpha,               out = '*';
    else,                           out = 'ns';
    end
end