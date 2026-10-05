function [sigPairs, statsTbl, anovaTbl, workTbl] = runOneWayANOVA(tbl, metricVec, options)
%RUNONEWAYANOVA  One-way ANOVA with optional bioRep averaging
%
% Fits standard one-way ANOVA:
%   metric ~ condition
%
% Optionally collapses technical replicates:
%   MeanBioRep = true → bioRep means used (recommended sanity check)
%
% -------------------------------------------------------------------------
% OPTIONS (Name-Value)
%   'CombinedVec'      : condition vector (e.g. 'CT_rest')
%   'BioRepCol'        : column for biological replicate
%   'ConditionLevels'  : ordered condition levels
%   'MeanBioRep'       : true/false (default false)
%   'Correction'       : 'fdr' | 'bonferroni' | 'none'
%   'Alpha'            : significance threshold (default 0.05)
%   'Verbose'          : print output
% -------------------------------------------------------------------------

arguments
    tbl table
    metricVec (:,1) double

    options.CombinedVec           = {}
    options.BioRepCol string      = "expCat"
    options.ConditionLevels cell  = {}
    options.MeanBioRep (1,1) logical = false
    options.Correction string     = "fdr"
    options.Alpha (1,1) double    = 0.05
    options.Verbose (1,1) logical = true
end

alphaVal = options.Alpha;

%% --- Build condition + bioRep vectors ---
bioVec  = string(tbl.(options.BioRepCol));
condVec = string(options.CombinedVec);

if isempty(options.ConditionLevels)
    condLevels = cellstr(unique(condVec(~isnan(metricVec)), 'stable'));
else
    condLevels = options.ConditionLevels;
end

%% --- Build working table ---
validMask = ~isnan(metricVec) & ismember(condVec, condLevels);

workTbl = table();
workTbl.metric    = metricVec(validMask);
workTbl.condition = categorical(condVec(validMask), condLevels);
workTbl.bioRep    = categorical(bioVec(validMask));

%% --- Optional: collapse to bioRep means ---
if options.MeanBioRep
    bioLevels  = categories(workTbl.bioRep);
    condLevelsPresent = categories(workTbl.condition);

    rows = {};
    for b = 1:numel(bioLevels)
        for k = 1:numel(condLevelsPresent)
            mask = workTbl.bioRep == bioLevels{b} & ...
                   workTbl.condition == condLevelsPresent{k};

            vals = workTbl.metric(mask);
            if isempty(vals), continue; end

            rows{end+1} = {mean(vals,'omitnan'), condLevelsPresent{k}, bioLevels{b}}; %#ok<AGROW>
        end
    end

    if isempty(rows)
        error('runOneWayANOVA: no data left after aggregation.');
    end
    tmp = vertcat(rows{:});
    % REBUILD TABLE instead of overwriting columns
    workTbl = table();
    workTbl.metric    = cell2mat(tmp(:,1));
    workTbl.condition = categorical(tmp(:,2), condLevels);
    workTbl.bioRep    = categorical(tmp(:,3));

    if options.Verbose
        fprintf('After bioRep averaging: %d rows\n', height(workTbl));
    end
end

%% --- Run one-way ANOVA ---
[p, anovaTbl, stats] = anova1(workTbl.metric, workTbl.condition, 'off');

if options.Verbose
    fprintf('\n=== One-way ANOVA ===\n');
    fprintf('p = %.4g\n', p);
end

%% --- Pairwise comparisons (Tukey) ---
[c,~,~,gnames] = multcompare(stats, 'Display', 'off', ...
    'CriticalValueType', 'lsd');  % LSD = uncorrected pairwise t-tests
rawP = c(:,6);

% Then your correction switch block is valid
nPairs = size(c,1);
comboA = gnames(c(:,1));
comboB = gnames(c(:,2));

% rawP = c(:,6);

%% --- Multiple comparison correction ---
switch lower(options.Correction)
    case 'fdr'
        adjP = mafdr(rawP, 'BHFDR', true);
    case 'bonferroni'
        adjP = min(rawP * nPairs, 1);
    otherwise
        adjP = rawP;
end

%% --- Build outputs ---
sigPairs = {};
sigLabel = cell(nPairs,1);

for i = 1:nPairs
    sigLabel{i} = pToStars(adjP(i), alphaVal);
    if adjP(i) < alphaVal
        sigPairs{end+1} = {comboA{i}, comboB{i}, sigLabel{i}}; %#ok<AGROW>
    end
end

statsTbl = table(comboA, comboB, rawP, adjP, sigLabel, ...
    'VariableNames', {'ComboA','ComboB','pRaw','pAdj','Sig'});

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

function out = pToStars(p, alpha)
    if nargin < 2, alpha = 0.05; end
    if     p < alpha * 0.001/0.05, out = '***';
    elseif p < alpha * 0.01/0.05,  out = '**';
    elseif p < alpha,               out = '*';
    else,                           out = 'ns';
    end
end