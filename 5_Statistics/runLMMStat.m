function [sigPairs, statsTbl, lmeModel, workTbl, dodTbl] = runLMMStat(tbl, metricVec, options)
%RUNLMMSTAT  Linear Mixed Model (LMM) statistical analysis for grouped,
%            repeated-measures data with unbalanced designs.
%
%   Fits:  metric ~ condition + (1|bioRep)
%
%   'condition' encodes all levels directly (e.g. CT_rest, CCM_1kPa),
%   so no separate group factor is needed. bioRep random intercept
%   accounts for unbalanced biological replicates across conditions.
%
% -------------------------------------------------------------------------
% USAGE
%   [sigPairs, statsTbl, lmeModel, workTbl] = runLLMStat(tbl, metricVec)
%   [sigPairs, statsTbl, lmeModel, workTbl] = runLLMStat(tbl, metricVec, Name, Value, ...)
%
% REQUIRED INPUTS
%   tbl        : MATLAB table containing all relevant columns.
%   metricVec  : Numeric vector whose rows correspond 1-to-1 with tbl.
%
% OPTIONAL NAME-VALUE INPUTS
%   'ConditionCol'   : Column name for the condition factor (all levels).
%                      Default: "extrainfoCat"
%   'BioRepCol'      : Column name for biological replicates.
%                      Default: "expCat"
%   'ConditionLevels': Cell array of strings — ordered levels for condition.
%                      Default: unique values found in ConditionCol.
%   'Contrasts'      : 'all' (default) | cell array of {keyA,keyB} pairs.
%   'Correction'     : 'fdr' (default) | 'bonferroni' | 'none'
%   'Alpha'          : Significance threshold (post-correction). Default 0.05.
%   'Verbose'        : Print LMM summary + contrast table. Default: true.cellMeanVecs
%   'DiffOfDiffs'    : {} (default, off) | 'auto' | explicit cell array.
%                      Difference-of-differences (interaction) contrasts,
%                      tested in the SAME fitted model (no refit):
%                         (G2_B - G2_A) - (G1_B - G1_A)
%                      e.g. (CCM_10kPa - CCM_1kPa) - (CT_10kPa - CT_1kPa).
%                      'auto' builds them from labels 'GENOTYPE_LEVEL_TYPE'.
%                      Explicit form: { {'CCM_10kPa_TS','CCM_1kPa_TS', ...
%                                        'CT_10kPa_TS','CT_1kPa_TS'}, ... }
%                      (order = G2_B, G2_A, G1_B, G1_A)
%   'DiffOfDiffsGroups'          : {reference, other} genotype. Default {'CT','CCM'}
%   'DiffOfDiffsLevelOrder'      : level order; later minus earlier.
%                                  Default {'rest','R1','R2','1kPa','10kPa'}
%   'DiffOfDiffsExcludePatterns' : level pairs to skip in 'auto'.
%                                  Default {{'rest','kPa'}} (glass vs gel)
%                      ContrastExcludePatterns/ContrastRemove are also respected.
%   'DiffOfDiffsCorrection'      : 'separate' (default) | 'joint'
%                      'separate': the DoD contrasts form their own FDR family,
%                                  so all existing pAdj/sigPairs are unchanged.
%                      'joint'   : DoD and pairwise contrasts share one family
%                                  (existing pAdj/sigPairs may change).
%
% OUTPUTS
%   sigPairs   : { {keyA,keyB,'***'}, {keyA,keyC,'*'}, ... } — significant pairs.
%   statsTbl   : Table with one row per contrast (see below).
%   lmeModel   : Fitted LinearMixedModel object (for diagnostics).
%   workTbl    : Filtered working table passed to fitlme.
%   dodTbl     : Difference-of-differences results (empty table if not requested).
%                DoD contrasts are NOT added to sigPairs (brackets span 2 bars).

% =========================================================================
% 1. Arguments
% =========================================================================
arguments
    tbl                         table
    metricVec                   (:,1) double

    options.CombinedVec                 = {}
    options.ConditionCol        string  = "extrainfoCat"
    options.BioRepCol           string  = "expCat"
    options.ConditionLevels     cell    = {}
    options.Contrasts                   = 'all'
    options.ContrastRemove      cell    = {}
    options.ContrastRemoveG1    cell    = {}
    options.ContrastRemoveG2    cell    = {}
    options.ContrastExcludePatterns   cell    = {} % e.g. {{'kPa','R1'}, {'kPa','R2'}}
    options.DiffOfDiffs                 = {}
    options.DiffOfDiffsGroups   cell    = {'CT','CCM'}
    options.DiffOfDiffsLevelOrder cell  = {'rest','R1','R2','1kPa','10kPa'}
    options.DiffOfDiffsExcludePatterns cell = {{'rest','kPa'}}
    options.DiffOfDiffsCorrection string {mustBeMember(options.DiffOfDiffsCorrection, ["separate","joint"])} = "separate"
    options.Correction          string  = "fdr"
    options.CollapseRep         (1,1) logical = true
    options.Alpha               (1,1) double ...
        {mustBeGreaterThan(options.Alpha,0), mustBeLessThan(options.Alpha,1)} = 0.05
    options.Verbose             (1,1) logical = true
end

condCol  = char(options.ConditionCol);
bioCol   = char(options.BioRepCol);
alphaVal = options.Alpha;
dodTbl   = table();

% =========================================================================
% 2. Validate inputs
% =========================================================================
assert(height(tbl) == numel(metricVec), ...
    'runLLMStat: tbl height (%d) must match numel(metricVec) (%d).', ...
    height(tbl), numel(metricVec));

assert(ismember(bioCol, tbl.Properties.VariableNames), ...
    'runLLMStat: column "%s" not found in tbl.', bioCol);

if isempty(options.CombinedVec)
    assert(ismember(condCol, tbl.Properties.VariableNames), ...
        'runLLMStat: column "%s" not found in tbl.', condCol);
end

% =========================================================================
% 3. Resolve condition levels
% =========================================================================
bioVec = string(tbl.(char(options.BioRepCol)));

if ~isempty(options.CombinedVec)
    condVec = string(options.CombinedVec);
else
    condVec = string(tbl.(char(options.ConditionCol)));
end

if isempty(options.ConditionLevels)
    condLevels = cellstr(unique(condVec(~isnan(metricVec) & condVec ~= ""), 'stable'));
else
    condLevels = options.ConditionLevels;
end

nCond = numel(condLevels);

% =========================================================================
% 4. Build working table (valid rows only)
% =========================================================================
validMask = ~isnan(metricVec) & ismember(condVec, condLevels);

workTbl           = table();
workTbl.metric    = metricVec(validMask);
workTbl.condition = categorical(condVec(validMask), condLevels, 'Ordinal', false);
workTbl.bioRep    = categorical(bioVec(validMask));

if height(workTbl) == 0
    warning('runLLMStat: no valid rows after filtering.');
    [sigPairs, statsTbl, lmeModel] = deal({}, table(), []);
    return
end

if nCond < 2
    warning('runLLMStat: fewer than 2 condition levels; nothing to test.');
    [sigPairs, statsTbl, lmeModel] = deal({}, table(), []);
    return
end

if options.CollapseRep
    % 4b. Aggregate to bioRep means
    %     Collapses technical replicates — the bioRep is now the unit of analysis.
    %     This matches the statistical intent of the SuperPlots approach but
    %     with proper LMM contrasts and FDR correction on top.
    bioRepLevels = categories(workTbl.bioRep);
    condLevelsPresent = categories(workTbl.condition);
    
    aggRows = {};
    for b = 1:numel(bioRepLevels)
        for k = 1:numel(condLevelsPresent)
            mask = (string(workTbl.bioRep)    == bioRepLevels{b}) & ...
                   (string(workTbl.condition) == condLevelsPresent{k});
            vals = workTbl.metric(mask);
            if sum(~isnan(vals)) == 0, continue; end   % this bioRep absent from this condition
    
            aggRows{end+1} = {mean(vals,'omitnan'), condLevelsPresent{k}, bioRepLevels{b}}; %#ok<AGROW>
        end
    end
    
    % Rebuild workTbl from aggregated rows
    aggMat        = vertcat(aggRows{:});
    workTbl       = table();
    workTbl.metric    = cell2mat(aggMat(:,1));
    workTbl.condition = categorical(aggMat(:,2), condLevels, 'Ordinal', false);
    workTbl.bioRep    = categorical(aggMat(:,3));
    
    if options.Verbose
        fprintf('  After aggregation to bioRep means: %d rows\n', height(workTbl));
    end
end

% =========================================================================
% 5. Fit LMM:  metric ~ condition + (1|bioRep)
% =========================================================================
if options.Verbose
    fprintf('\n=== runLLMStat: fitting LMM ===\n');
    fprintf('  metric ~ condition + (1|bioRep)\n');
    fprintf('  Conditions (%d): %s\n', nCond, strjoin(condLevels, ', '));
    fprintf('  Valid obs      : %d\n', height(workTbl));
    fprintf('  Bio replicates : %d\n', numel(unique(workTbl.bioRep)));
end

try
    lmeModel = fitlme(workTbl, ...
        'metric ~ condition + (1|bioRep)', ...
        'FitMethod',       'REML', ...
        'DummyVarCoding',  'effects');
catch ME
    error('runLLMStat: fitlme failed — %s', ME.message);
end

if options.Verbose
    disp(lmeModel);
end

% =========================================================================
% 6. Enumerate contrasts
%    comboKeys == condLevels (no group prefix needed)
% =========================================================================
comboKeys = condLevels;   % e.g. {'CT_rest','CT_R1','CCM_rest', ...}

if ischar(options.Contrasts) || isstring(options.Contrasts)
    switch lower(char(options.Contrasts))
        case 'all'
            pairsToTest = allPairs(comboKeys);
        otherwise
            error('runLLMStat: unknown Contrasts mode "%s". Use ''all'' or a cell array.', ...
                options.Contrasts);
    end
elseif iscell(options.Contrasts)
    pairsToTest = options.Contrasts;
    for r = 1:numel(pairsToTest)
        p = pairsToTest{r};
        assert(iscell(p) && numel(p) == 2, ...
            'runLLMStat: each Contrasts element must be a 1×2 cell {keyA, keyB}.');
        for side = 1:2
            assert(ismember(p{side}, comboKeys), ...
                'runLLMStat: key "%s" not in condLevels.\nAvailable: %s', ...
                p{side}, strjoin(comboKeys, ', '));
        end
    end
else
    error('runLLMStat: Contrasts must be ''all'' or a cell array of pairs.');
end

% =========================================================================
% 6b. Build unified exclusion set and filter pairsToTest
% =========================================================================

% % ── Remove combinations of conditions
% if ~isempty(options.ContrastRemoveG1) && ~isempty(options.ContrastRemoveG2)
%     options.ContrastRemove = {};
%     for i = 1:numel(options.ContrastRemoveG1)
%         for j = 1:numel(options.ContrastRemoveG2)
%             options.ContrastRemove{end+1} = {options.ContrastRemoveG1{i}, options.ContrastRemoveG2{j}};
%         end
%     end
% end
% 
% % ── Remove explicitly excluded pairs
% if ~isempty(options.ContrastRemove)
%     removeSet = options.ContrastRemove;
%     keepMask  = true(numel(pairsToTest), 1);
%     for ri = 1:numel(pairsToTest)
%         pA = pairsToTest{ri}{1};
%         pB = pairsToTest{ri}{2};
%         for rj = 1:numel(removeSet)
%             rA = removeSet{rj}{1};
%             rB = removeSet{rj}{2};
%             % Match in either order
%             if (strcmp(pA,rA) && strcmp(pB,rB)) || ...
%                (strcmp(pA,rB) && strcmp(pB,rA))
%                 keepMask(ri) = false;
%                 break
%             end
%         end
%     end
%     pairsToTest = pairsToTest(keepMask);
% end

excludePairs = {};  % will collect all {keyA, keyB} pairs to remove

% ── From ContrastRemoveG1 × ContrastRemoveG2 cross-product
if ~isempty(options.ContrastRemoveG1) && ~isempty(options.ContrastRemoveG2)
    for i = 1:numel(options.ContrastRemoveG1)
        for j = 1:numel(options.ContrastRemoveG2)
            excludePairs{end+1} = {options.ContrastRemoveG1{i}, options.ContrastRemoveG2{j}}; %#ok<AGROW>
        end
    end
end

% ── From explicit ContrastRemove list
for i = 1:numel(options.ContrastRemove)
    excludePairs{end+1} = options.ContrastRemove{i}; %#ok<AGROW>
end

% ── From pattern-based exclusions: expand against comboKeys
if ~isempty(options.ContrastExcludePatterns)
    for pi = 1:numel(options.ContrastExcludePatterns)
        patA = options.ContrastExcludePatterns{pi}{1};
        patB = options.ContrastExcludePatterns{pi}{2};
        for i = 1:numel(comboKeys)
            for j = 1:numel(comboKeys)
                if i == j, continue; end
                if contains(comboKeys{i}, patA) && contains(comboKeys{j}, patB)
                    excludePairs{end+1} = {comboKeys{i}, comboKeys{j}}; %#ok<AGROW>
                end
            end
        end
    end
end

% ── Canonicalise: sort each pair alphabetically, then deduplicate
canonExclude = cellfun(@(p) sort(p), excludePairs, 'UniformOutput', false);
canonExclude = unique(cellfun(@(p) strjoin(p, '|||'), canonExclude, 'UniformOutput', false));

% ── Apply to pairsToTest
keepMask = true(numel(pairsToTest), 1);
for ri = 1:numel(pairsToTest)
    canon = strjoin(sort(pairsToTest{ri}), '|||');
    if ismember(canon, canonExclude)
        keepMask(ri) = false;
    end
end
pairsToTest = pairsToTest(keepMask);

% ── Check pairsToTest
nPairs = numel(pairsToTest);
if nPairs == 0
    warning('runLLMStat: no contrasts to test.');
    [sigPairs, statsTbl] = deal({}, table());
    return
end

if options.Verbose
    fprintf('  Contrasts to test: %d\n', nPairs);
end

% =========================================================================
% 7. Per-cell descriptive statistics
% =========================================================================
cellStats = struct('key',{}, 'meanVal',{}, 'nObs',{}, 'nBioRep',{});
for k = 1:nCond
    mask = string(workTbl.condition) == condLevels{k};
    vals = workTbl.metric(mask);
    bios = workTbl.bioRep(mask);

    cellStats(k).key     = condLevels{k};
    cellStats(k).meanVal = mean(vals, 'omitnan');
    cellStats(k).nObs    = sum(~isnan(vals));
    cellStats(k).nBioRep = numel(unique(bios(~isnan(vals))));
end

% =========================================================================
% 8. Contrast vectors via fixed-effects design matrix
%
%    For each condition level, take the mean of its design-matrix rows.
%    This gives a stable cell-mean vector even with effects coding,
%    and is robust to small within-cell imbalance in tech replicates.
% =========================================================================
X    = designMatrix(lmeModel, 'Fixed');   % N × nFE
nFE  = size(X, 2);

cellMeanVecs = zeros(nCond, nFE);
for k = 1:nCond
    rowMask = string(workTbl.condition) == condLevels{k};
    if any(rowMask)
        % Average the design-matrix rows for this condition —
        % more robust than picking a single representative row
        cellMeanVecs(k, :) = mean(X(rowMask, :), 1);
    else
        warning('runLLMStat: no rows for condition "%s".', condLevels{k});
    end
end

betaHat  = fixedEffects(lmeModel);
CovBeta  = lmeModel.CoefficientCovariance;

rawP      = nan(nPairs, 1);
estimates = nan(nPairs, 1);
seVals    = nan(nPairs, 1);
tStats    = nan(nPairs, 1);
dfVals    = nan(nPairs, 1);

for r = 1:nPairs
    keyA = pairsToTest{r}{1};
    keyB = pairsToTest{r}{2};
    idxA = find(strcmp(comboKeys, keyA));
    idxB = find(strcmp(comboKeys, keyB));

    L = cellMeanVecs(idxA, :) - cellMeanVecs(idxB, :);

    [pVal, ~, ~]  = coefTest(lmeModel, L);
    rawP(r)       = pVal;
    estimates(r)  = L * betaHat;
    seVals(r)     = sqrt(max(0, L * CovBeta * L'));
    tStats(r)     = estimates(r) / seVals(r);
    dfVals(r)     = lmeModel.DFE;
end

% =========================================================================
% 9. Multiple-comparison correction
% =========================================================================
switch lower(char(options.Correction))
    case 'fdr'
        adjP = mafdr(rawP, 'BHFDR', true);
    case 'bonferroni'
        adjP = min(rawP * nPairs, 1);
    case 'none'
        adjP = rawP;
    otherwise
        warning('runLLMStat: unknown correction "%s"; using FDR.', options.Correction);
        adjP = mafdr(rawP, 'BHFDR', true);
end

% -------------------------------------------------------------------------
% 9b. Difference-of-differences contrasts (same model, no refit)
% -------------------------------------------------------------------------
dodDefs = {};   % each: {G2_B, G2_A, G1_B, G1_A, label}
dd = options.DiffOfDiffs;
if iscell(dd) && ~isempty(dd)
    for q = 1:numel(dd)
        d = dd{q};
        assert(iscell(d) && numel(d) == 4 && all(ismember(d, comboKeys)), ...
            ['runLLMStat: each DiffOfDiffs entry must be a 1x4 cell of keys present ' ...
             'in condLevels: {G2_B, G2_A, G1_B, G1_A}.']);
        dodDefs{end+1} = [d(:)', {sprintf('(%s - %s) - (%s - %s)', d{1}, d{2}, d{3}, d{4})}]; %#ok<AGROW>
    end
elseif (ischar(dd) || isstring(dd)) && strcmpi(dd, 'auto')
    dodDefs = autoDoD(comboKeys, options.DiffOfDiffsGroups, ...
        options.DiffOfDiffsLevelOrder, canonExclude, options.DiffOfDiffsExcludePatterns);
end

% keep only contrasts whose four cells contain data
if ~isempty(dodDefs)
    hasAll = cellfun(@(d) all(arrayfun(@(k) ...
        cellStats(strcmp({cellStats.key}, d{k})).nObs > 0, 1:4)), dodDefs);
    dodDefs = dodDefs(hasAll);
end
nDoD = numel(dodDefs);

if nDoD > 0
    dodEst = nan(nDoD,1); dodSE = nan(nDoD,1); dodT = nan(nDoD,1);
    dodDF  = nan(nDoD,1); dodRawP = nan(nDoD,1);
    for q = 1:nDoD
        d  = dodDefs{q};
        i1 = find(strcmp(comboKeys, d{1}));   % G2_B
        i2 = find(strcmp(comboKeys, d{2}));   % G2_A
        i3 = find(strcmp(comboKeys, d{3}));   % G1_B
        i4 = find(strcmp(comboKeys, d{4}));   % G1_A
        L  = (cellMeanVecs(i1,:) - cellMeanVecs(i2,:)) - ...
             (cellMeanVecs(i3,:) - cellMeanVecs(i4,:));
        [pVal, ~, ~] = coefTest(lmeModel, L);
        dodRawP(q) = pVal;
        dodEst(q)  = L * betaHat;
        dodSE(q)   = sqrt(max(0, L * CovBeta * L'));
        dodT(q)    = dodEst(q) / dodSE(q);
        dodDF(q)   = lmeModel.DFE;
    end

    switch lower(char(options.DiffOfDiffsCorrection))
        case 'separate'
            dodAdjP = adjustP(dodRawP, options.Correction);
        case 'joint'
            allAdj  = adjustP([rawP; dodRawP], options.Correction);
            adjP    = allAdj(1:nPairs);          % existing contrasts re-adjusted
            dodAdjP = allAdj(nPairs+1:end);
    end

    dodTbl = table(string(cellfun(@(d) d{5}, dodDefs(:), 'UniformOutput', false)), ...
        dodEst, dodSE, dodT, dodDF, dodRawP, dodAdjP, ...
        cellfun(@pToStars, num2cell(dodAdjP), 'UniformOutput', false), ...
        'VariableNames', {'Contrast','Estimate','SE','tStat','DF','pRaw','pAdj','Sig'});
end

% =========================================================================
% 10. Build statsTbl & sigPairs
% =========================================================================
comboA   = cell(nPairs,1);  comboB   = cell(nPairs,1);
meanA    = nan(nPairs,1);   meanB    = nan(nPairs,1);
nObsA    = nan(nPairs,1);   nObsB    = nan(nPairs,1);
nBioA    = nan(nPairs,1);   nBioB    = nan(nPairs,1);
sigLabel = cell(nPairs,1);

for r = 1:nPairs
    keyA = pairsToTest{r}{1};
    keyB = pairsToTest{r}{2};
    comboA{r} = keyA;
    comboB{r} = keyB;

    sA = cellStats(strcmp({cellStats.key}, keyA));
    sB = cellStats(strcmp({cellStats.key}, keyB));
    meanA(r)    = sA.meanVal;    meanB(r)    = sB.meanVal;
    nObsA(r)    = sA.nObs;       nObsB(r)    = sB.nObs;
    nBioA(r)    = sA.nBioRep;    nBioB(r)    = sB.nBioRep;
    sigLabel{r} = pToStars(adjP(r)); %, alphaVal);
end

statsTbl = table(comboA, comboB, estimates, seVals, tStats, dfVals, ...
                 rawP, adjP, sigLabel, ...
                 meanA, meanB, nObsA, nObsB, nBioA, nBioB, ...
    'VariableNames', {'ComboA','ComboB','Estimate','SE','tStat','DF', ...
                      'pRaw','pAdj','Sig', ...
                      'MeanA','MeanB','NObsA','NObsB','NBioRepA','NBioRepB'});

sigPairs = {};
for r = 1:nPairs
    if adjP(r) < alphaVal
        sigPairs{end+1} = {comboA{r}, comboB{r}, sigLabel{r}}; %#ok<AGROW>
    end
end

% =========================================================================
% 11. Verbose output
% =========================================================================
if options.Verbose
    fprintf('\n--- Contrast results (%s correction, α=%.2f) ---\n', ...
        options.Correction, alphaVal);
    fprintf('%-22s  %-22s  %8s  %8s  %6s\n', 'ComboA','ComboB','pRaw','pAdj','Sig');
    fprintf('%s\n', repmat('-',1,72));
    for r = 1:nPairs
        fprintf('%-22s  %-22s  %8.4g  %8.4g  %s\n', ...
            comboA{r}, comboB{r}, rawP(r), adjP(r), sigLabel{r});
    end
    fprintf('\n  Significant pairs: %d / %d\n\n', numel(sigPairs), nPairs);
    if nDoD > 0
        fprintf('--- Difference-of-differences contrasts (%s FDR family) ---\n', ...
            options.DiffOfDiffsCorrection);
        fprintf('%-62s  %9s  %8s  %8s  %s\n', 'Contrast','Estimate','pRaw','pAdj','Sig');
        fprintf('%s\n', repmat('-',1,98));
        for r = 1:nDoD
            fprintf('%-62s  %9.4g  %8.4g  %8.4g  %s\n', dodTbl.Contrast(r), ...
                dodTbl.Estimate(r), dodTbl.pRaw(r), dodTbl.pAdj(r), dodTbl.Sig{r});
        end
        fprintf('\n');
    end
end

end   % ── end main ──────────────────────────────────────────────────────


% =========================================================================
% LOCAL HELPERS
% =========================================================================

function pairs = allPairs(keys)
    n = numel(keys);
    pairs = {};
    for i = 1:n
        for j = i+1:n
            pairs{end+1} = {keys{i}, keys{j}}; %#ok<AGROW>
        end
    end
end

% function out = pToStars(p, alpha)
%     if nargin < 2, alpha = 0.05; end
%     if     p < alpha * 0.001/0.05, out = '***';
%     elseif p < alpha * 0.01/0.05,  out = '**';
%     elseif p < alpha,               out = '*';
%     else,                           out = 'ns';
%     end
% end
function out = pToStars(p)
    if p < 0.001,     out = '***';
    elseif p < 0.01,  out = '**';
    elseif p < 0.05,  out = '*';
    else,             out = 'ns';
    end
end

function adj = adjustP(p, method)
% Multiple-comparison correction, same options as the main function
    switch lower(char(method))
        case 'fdr',        adj = mafdr(p, 'BHFDR', true);
        case 'bonferroni', adj = min(p * numel(p), 1);
        case 'none',       adj = p;
        otherwise
            warning('runLLMStat: unknown correction "%s"; using FDR.', method);
            adj = mafdr(p, 'BHFDR', true);
    end
end

function defs = autoDoD(keys, grp, lvlOrder, canonExclude, excl)
% Build (G2_B - G2_A) - (G1_B - G1_A) for every level pair (A earlier than B in
% lvlOrder) where all four cells exist. Keys look like 'CCM_10kPa_TS'.
    defs = {};
    g1 = grp{1};  g2 = grp{2};
    canon = @(x, y) strjoin(sort({x, y}), '|||');
    for a = 1:numel(lvlOrder)
        for b = a+1:numel(lvlOrder)
            lvA = lvlOrder{a};  lvB = lvlOrder{b};
            skip = false;
            for e = 1:numel(excl)
                pA = excl{e}{1};  pB = excl{e}{2};
                if (contains(lvA,pA) && contains(lvB,pB)) || ...
                   (contains(lvA,pB) && contains(lvB,pA))
                    skip = true;
                end
            end
            if skip, continue; end

            prefix = [g1 '_' lvA];
            for k = 1:numel(keys)
                if ~startsWith(keys{k}, prefix), continue; end
                rest = keys{k}(numel(prefix)+1:end);        % e.g. '_TS'
                if ~isempty(rest) && rest(1) ~= '_', continue; end
                k11A = [g1 '_' lvA rest];  k11B = [g1 '_' lvB rest];
                k21A = [g2 '_' lvA rest];  k21B = [g2 '_' lvB rest];
                if ~all(ismember({k11A,k11B,k21A,k21B}, keys)), continue; end
                % respect exclusions applied to the underlying pairwise contrasts
                if ismember(canon(k11A,k11B), canonExclude) || ...
                   ismember(canon(k21A,k21B), canonExclude)
                    continue
                end
                defs{end+1} = {k21B, k21A, k11B, k11A, ...
                    sprintf('(%s - %s) - (%s - %s)', k21B, k21A, k11B, k11A)}; %#ok<AGROW>
            end
        end
    end
end