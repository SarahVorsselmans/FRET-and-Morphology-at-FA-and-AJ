function [sigPairs, statsTbl, anovaTbl, workTbl] = runNestedANOVAStat(tbl, metricVec, options)
%RUNNESTEDANOVASTAT  Nested ANOVA statistical analysis for grouped,
%                    repeated-measures data where biological replicates are
%                    nested within condition × group cells.
%
%   Model hierarchy (Type III SS, balanced or unbalanced via anovan):
%
%       Level 1  (fixed)  : condition
%       Level 1  (fixed)  : group
%       Level 1  (fixed)  : condition × group  interaction
%       Level 2  (random) : bioRep(condition × group)   ← nested error term
%       Level 3  (error)  : residual (technical replicates within bioRep)
%
%   Fixed effects are tested against the nested bioRep MS (not the residual),
%   which correctly accounts for pseudo-replication when multiple technical
%   measurements come from the same biological sample.
%
%   Post-hoc pairwise contrasts between condition×group cells use a
%   Tukey-Kramer-style t-test with:
%       SE  = sqrt( MS_bioRep_nested * (1/n_A + 1/n_B) )
%       DF  = DF_bioRep_nested   (conservative Satterthwaite not needed
%                                 because the same error term is used for
%                                 all contrasts)
%
%   This differs from the LMM approach (runLMMStat) in the following ways:
%     • No random-effect parameters are estimated via REML; the nested
%       variance is captured implicitly through the MS ratio.
%     • Missing biological replicates across conditions cause listwise
%       exclusion of those bio-rep IDs (like rmANOVA), so unbalanced
%       designs with many missing cells are better handled by the LMM.
%     • The ANOVA table is explicit and directly readable, which is often
%       required by reviewers in biology journals.
%
% -------------------------------------------------------------------------
% USAGE
%   [sigPairs, statsTbl, anovaTbl] = runNestedANOVAStat(tbl, metricVec)
%   [sigPairs, statsTbl, anovaTbl] = runNestedANOVAStat(tbl, metricVec, Name, Value, ...)
%
% REQUIRED INPUTS
%   tbl        : MATLAB table containing all relevant columns.
%   metricVec  : Numeric vector whose rows correspond 1-to-1 with tbl.
%
% OPTIONAL NAME-VALUE INPUTS
%   'ConditionCol'   : Column name for the condition factor.
%                      Default: "extrainfoCat"
%   'GroupCol'       : Column name for the group factor (e.g. CT/CCM).
%                      Default: "diseaseCat"
%   'BioRepCol'      : Column name for biological replicates.
%                      Default: "expCat"
%   'ConditionLevels': Ordered levels for condition (cell array of strings).
%                      Default: unique values found in ConditionCol.
%   'GroupLevels'    : Which group values to include (cell array of strings).
%                      Default: all non-empty unique values in GroupCol.
%
%   'Contrasts'      : What pairs to test. One of:
%                        'all'       – all pairwise condition×group combos (default)
%                        'within'    – only group comparisons within each condition
%                        cell array  – explicit {{keyA,keyB}, {keyC,keyD}, ...}
%
%   'Correction'     : Multiple-comparison correction.
%                        'fdr'        – Benjamini-Hochberg (default)
%                        'bonferroni'
%                        'none'
%   'Alpha'          : Significance threshold. Default: 0.05.
%   'Verbose'        : Print ANOVA table and contrast results. Default: true.
%   'KeySep'         : Separator for combo-keys. Default: "_".
%
% OUTPUTS
%   sigPairs   : Cell array of significant pairs compatible with
%                plotMetricScatter bracket-drawing code:
%                  { {keyA, keyB, '***'}, ... }
%
%   statsTbl   : Table with one row per tested contrast:
%                  ComboA, ComboB, Estimate, SE, tStat, DF, pRaw, pAdj,
%                  Sig, MeanA, MeanB, NObsA, NObsB, NBioRepA, NBioRepB
%
%   anovaTbl   : ANOVA table returned by anovan (useful for reporting and
%                checking that the nested bioRep term is significant,
%                which would validate the nested design choice).
%
%   workTbl    : Filtered working table used for the analysis.
%
% DEPENDENCIES
%   Statistics and Machine Learning Toolbox (anovan, multcompare)
%
% NOTES
%   anovan is called with:
%       'random'  = index of bioRep factor  (makes it a random factor)
%       'nested'  = nesting matrix          (bioRep nested in cond×group)
%       'varnames'= {'condition','group','bioRep'}
%       'model'   = 'full' restricted to the three terms above plus interaction
%   The MS for the nested bioRep term is extracted from the anovaTbl and
%   used as the denominator for all fixed-effect F-tests and post-hoc SEs.

% =========================================================================
% 1. Arguments
% =========================================================================
arguments
    tbl                         table
    metricVec                   (:,1) double

    options.ConditionCol        string = "extrainfoCat"
    options.GroupCol            string = "diseaseCat"
    options.BioRepCol           string = "expCat"

    options.ConditionLevels     cell   = {}
    options.GroupLevels         cell   = {}

    options.Contrasts                  = 'all'

    options.Correction          string = "fdr"
    options.Alpha               (1,1) double ...
        {mustBeGreaterThan(options.Alpha, 0), mustBeLessThan(options.Alpha, 1)} = 0.05
    options.Verbose             (1,1) logical = true
    options.KeySep              string = "_"
end

condCol  = char(options.ConditionCol);
grpCol   = char(options.GroupCol);
bioCol   = char(options.BioRepCol);
sep      = char(options.KeySep);
alphaVal = options.Alpha;

% =========================================================================
% 2. Validate inputs
% =========================================================================
assert(height(tbl) == numel(metricVec), ...
    'runNestedANOVAStat: tbl height (%d) must match numel(metricVec) (%d).', ...
    height(tbl), numel(metricVec));

for c = {condCol, grpCol, bioCol}
    assert(ismember(c{1}, tbl.Properties.VariableNames), ...
        'runNestedANOVAStat: column "%s" not found in tbl.', c{1});
end

% =========================================================================
% 3. Resolve levels
% =========================================================================
condVec = string(tbl.(condCol));
grpVec  = string(tbl.(grpCol));
bioVec  = string(tbl.(bioCol));

if isempty(options.ConditionLevels)
    condLevels = cellstr(unique(condVec(~isnan(metricVec) & condVec ~= ""), 'stable'));
else
    condLevels = options.ConditionLevels;
end

if isempty(options.GroupLevels)
    grpLevels = cellstr(unique(grpVec(~isnan(metricVec) & grpVec ~= ""), 'stable'));
else
    grpLevels = options.GroupLevels;
end

% =========================================================================
% 4. Build working table (valid rows only)
% =========================================================================
validMask = ~isnan(metricVec) ...
    & ismember(condVec, condLevels) ...
    & ismember(grpVec,  grpLevels);

workTbl = table();
workTbl.metric    = metricVec(validMask);
workTbl.condition = categorical(condVec(validMask), condLevels, 'Ordinal', false);
workTbl.group     = categorical(grpVec(validMask),  grpLevels,  'Ordinal', false);
workTbl.bioRep    = categorical(bioVec(validMask));

% Create a combined cell ID: each unique condition×group×bioRep triplet
% gets a unique label.  This is what "nests" inside each cell.
workTbl.cellGrp = categorical( ...
    string(workTbl.group) + sep + string(workTbl.condition));

if height(workTbl) == 0
    warning('runNestedANOVAStat: no valid rows after filtering.');
    sigPairs = {}; statsTbl = table(); anovaTbl = table(); return;
end

nCond = numel(condLevels);
nGrp  = numel(grpLevels);
nCell = nCond * nGrp;

if nCell < 2
    warning('runNestedANOVAStat: fewer than 2 condition×group cells; nothing to test.');
    sigPairs = {}; statsTbl = table(); anovaTbl = table(); return;
end

% =========================================================================
% 5. Run nested ANOVA via anovan
%
%    Factor layout passed to anovan:
%      col 1 : condition  (fixed)
%      col 2 : group      (fixed)
%      col 3 : bioRep     (random, nested within condition × group)
%
%    Nesting matrix N (3×3):
%      N(i,j) = 1 means factor i is nested in factor j.
%      bioRep (row 3) is nested in condition (col 1) AND group (col 2).
%      No other nesting.
%
%    'model' string encodes which terms appear:
%      [1 0 0]  condition main effect
%      [0 1 0]  group main effect
%      [1 1 0]  condition × group interaction
%      [0 0 1]  bioRep main effect  (= bioRep(cond×group) because of nesting)
%    We deliberately omit higher-order terms with bioRep as they are not
%    estimable with typical biological replicate counts.
% =========================================================================
y        = workTbl.metric;
gCond    = cellstr(workTbl.condition);
gGrp     = cellstr(workTbl.group);
gBioRep  = cellstr(workTbl.bioRep);   % bioRep IDs are already unique per cell

% Nesting matrix: bioRep (factor 3) nested in condition (1) and group (2)
nestMat = [0 0 0;
           0 0 0;
           1 1 0];

modelMat = [1 0 0;   % condition
            0 1 0;   % group
            1 1 0;   % condition × group
            0 0 1];  % bioRep(condition × group)

if options.Verbose
    fprintf('\n=== runNestedANOVAStat: fitting nested ANOVA ===\n');
    fprintf('  Model: metric ~ condition * group + bioRep(condition*group)\n');
    fprintf('  Conditions (%d): %s\n', nCond, strjoin(condLevels, ', '));
    fprintf('  Groups     (%d): %s\n', nGrp,  strjoin(grpLevels,  ', '));
    fprintf('  Valid obs       : %d\n', height(workTbl));
    fprintf('  Bio replicates  : %d\n', numel(unique(workTbl.bioRep)));
end

try
    [pAnova, anovaTbl, anovaStats] = anovan( ...
        y, ...
        {gCond, gGrp, gBioRep}, ...
        'random',   3, ...
        'nested',   nestMat, ...
        'model',    modelMat, ...
        'varnames', {'condition', 'group', 'bioRep'}, ...
        'display',  'off');
catch ME
    error('runNestedANOVAStat: anovan failed — %s', ME.message);
end

if options.Verbose
    fprintf('\n--- ANOVA table ---\n');
    disp(anovaTbl);
end

% =========================================================================
% 6. Extract nested bioRep MS and DF for use in post-hoc tests
%
%    anovan returns anovaTbl as a cell array.  Row names are in column 1.
%    Locate the bioRep row (nested term) to get its MS and DF.
%    These serve as the error term for all pairwise t-tests below.
%
%    Column layout of anovaTbl (standard anovan output):
%      col 1 : Source name
%      col 2 : SS
%      col 3 : df
%      col 4 : MS
%      col 5 : F
%      col 6 : Prob>F
% =========================================================================
rowNames = anovaTbl(2:end-1, 1);   % skip header row and 'Total' row

header = anovaTbl(1,:);
msCol = find(strcmp(header, 'Mean Sq.'));
dfCol = find(strcmp(header, 'd.f.'));

% Find the bioRep nested term row (its name contains 'bioRep')
bioRepRowIdx = find(contains(string(rowNames), 'bioRep', 'IgnoreCase', true), 1);

if isempty(bioRepRowIdx)
    warning(['runNestedANOVAStat: could not locate bioRep nested term in ANOVA table. ' ...
             'Falling back to residual MS for post-hoc tests.']);
    % Fall back to Error row
    bioRepRowIdx = find(contains(string(rowNames), 'Error', 'IgnoreCase', true), 1);
end

% +1 offset because rowNames skips the header (anovaTbl row 1)
dataRow   = bioRepRowIdx + 1;
MS_nested = anovaTbl{dataRow, msCol};
DF_nested = anovaTbl{dataRow, dfCol};

if options.Verbose
    fprintf('  Nested bioRep MS = %.4g  (DF = %g)  — used as error term\n', ...
        MS_nested, DF_nested);
end

% =========================================================================
% 7. Build combo keys and enumerate contrasts
%    Key convention: "<group><sep><condition>"  e.g. "CT_rest"
% =========================================================================
comboKeys = cell(nCell, 1);
idx = 0;
for g = 1:nGrp
    for k = 1:nCond
        idx = idx + 1;
        comboKeys{idx} = [grpLevels{g}, sep, condLevels{k}];
    end
end

if ischar(options.Contrasts) || isstring(options.Contrasts)
    contrastMode = lower(char(options.Contrasts));
    switch contrastMode
        case 'all'
            pairsToTest = allPairs(comboKeys);

        case 'within'
            pairsToTest = {};
            for k = 1:nCond
                inCond = comboKeys(endsWith(comboKeys, [sep, condLevels{k}]));
                pairsToTest = [pairsToTest; allPairs(inCond)]; %#ok<AGROW>
            end

        otherwise
            error('runNestedANOVAStat: unknown Contrasts mode "%s".', contrastMode);
    end

elseif iscell(options.Contrasts)
    pairsToTest = options.Contrasts;
    for r = 1:numel(pairsToTest)
        p = pairsToTest{r};
        assert(iscell(p) && numel(p) == 2, ...
            'runNestedANOVAStat: each element of Contrasts must be a 1×2 cell.');
        for side = 1:2
            assert(ismember(p{side}, comboKeys), ...
                'runNestedANOVAStat: key "%s" not found. Available: %s', ...
                p{side}, strjoin(comboKeys, ', '));
        end
    end
else
    error('runNestedANOVAStat: Contrasts must be ''all'', ''within'', or a cell array.');
end

nPairs = numel(pairsToTest);
if nPairs == 0
    warning('runNestedANOVAStat: no contrasts to test.');
    sigPairs = {}; statsTbl = table(); return;
end

if options.Verbose
    fprintf('  Contrasts to test: %d\n', nPairs);
end

% =========================================================================
% 8. Per-cell descriptive statistics
% =========================================================================
cellStats = struct();
for ci = 1:nCell
    key   = comboKeys{ci};
    parts = strsplit(key, sep);
    grpStr  = parts{1};
    condStr = strjoin(parts(2:end), sep);

    mask = (string(workTbl.group) == grpStr) & ...
           (string(workTbl.condition) == condStr);
    vals = workTbl.metric(mask);
    bios = workTbl.bioRep(mask);

    cellStats(ci).key      = key;
    cellStats(ci).meanVal  = mean(vals, 'omitnan');
    cellStats(ci).nObs     = sum(~isnan(vals));
    cellStats(ci).nBioRep  = numel(unique(bios(~isnan(vals))));
end

% =========================================================================
% 9. Pairwise post-hoc t-tests using the nested bioRep MS as error
%
%    For each pair (A, B):
%      t = (meanA - meanB) / sqrt( MS_nested * (1/nA + 1/nB) )
%      p from t-distribution with DF_nested degrees of freedom
%
%    nA and nB here are observation counts (not bio-rep counts) because
%    MS_nested already reflects the between-bio-rep variance at the
%    observation scale when nested within cells.
%
%    If you prefer bio-rep-level means as the unit, set nA = nBioRepA etc.
%    and use bio-rep means as cellStats.meanVal (adjust above accordingly).
% =========================================================================
rawP      = nan(nPairs, 1);
estimates = nan(nPairs, 1);
seVals    = nan(nPairs, 1);
tStats    = nan(nPairs, 1);

for r = 1:nPairs
    keyA = pairsToTest{r}{1};
    keyB = pairsToTest{r}{2};

    sA = cellStats(strcmp({cellStats.key}, keyA));
    sB = cellStats(strcmp({cellStats.key}, keyB));

    nA = sA.nObs;
    nB = sB.nObs;

    if nA == 0 || nB == 0 || MS_nested <= 0
        rawP(r)      = NaN;
        estimates(r) = NaN;
        seVals(r)    = NaN;
        tStats(r)    = NaN;
        continue;
    end

    estimates(r) = sA.meanVal - sB.meanVal;
    seVals(r)    = sqrt(MS_nested * (1/nA + 1/nB));
    tStats(r)    = estimates(r) / seVals(r);

    % Two-tailed p from t-distribution with DF_nested
    rawP(r) = 2 * tcdf(-abs(tStats(r)), DF_nested);
end

% =========================================================================
% 10. Multiple-comparison correction
% =========================================================================
switch lower(char(options.Correction))
    case 'fdr'
        adjP = mafdr(rawP, 'BHFDR', true);
    case 'bonferroni'
        adjP = min(rawP * nPairs, 1);
    case 'none'
        adjP = rawP;
    otherwise
        warning('runNestedANOVAStat: unknown correction "%s"; using FDR.', options.Correction);
        adjP = mafdr(rawP, 'BHFDR', true);
end

% =========================================================================
% 11. Build statsTbl
% =========================================================================
comboA   = cell(nPairs, 1);
comboB   = cell(nPairs, 1);
meanA    = nan(nPairs, 1);
meanB    = nan(nPairs, 1);
nObsA    = nan(nPairs, 1);
nObsB    = nan(nPairs, 1);
nBioA    = nan(nPairs, 1);
nBioB    = nan(nPairs, 1);
sigLabel = cell(nPairs, 1);
dfVals   = repmat(DF_nested, nPairs, 1);

for r = 1:nPairs
    keyA = pairsToTest{r}{1};
    keyB = pairsToTest{r}{2};
    comboA{r} = keyA;
    comboB{r} = keyB;

    sA = cellStats(strcmp({cellStats.key}, keyA));
    sB = cellStats(strcmp({cellStats.key}, keyB));
    meanA(r)  = sA.meanVal;
    meanB(r)  = sB.meanVal;
    nObsA(r)  = sA.nObs;
    nObsB(r)  = sB.nObs;
    nBioA(r)  = sA.nBioRep;
    nBioB(r)  = sB.nBioRep;
    sigLabel{r} = pToStars(adjP(r), alphaVal);
end

statsTbl = table(comboA, comboB, estimates, seVals, tStats, dfVals, rawP, adjP, sigLabel, ...
    meanA, meanB, nObsA, nObsB, nBioA, nBioB, ...
    'VariableNames', {'ComboA','ComboB','Estimate','SE','tStat','DF', ...
                      'pRaw','pAdj','Sig', ...
                      'MeanA','MeanB','NObsA','NObsB','NBioRepA','NBioRepB'});

% =========================================================================
% 12. Build sigPairs output
% =========================================================================
sigPairs = {};
for r = 1:nPairs
    if ~isnan(adjP(r)) && adjP(r) < alphaVal
        sigPairs{end+1} = {comboA{r}, comboB{r}, sigLabel{r}}; %#ok<AGROW>
    end
end

% =========================================================================
% 13. Verbose output
% =========================================================================
if options.Verbose
    fprintf('\n--- Contrast results (%s correction, α=%.2f) ---\n', ...
        options.Correction, alphaVal);
    fprintf('%-22s  %-22s  %8s  %8s  %6s\n', 'ComboA', 'ComboB', 'pRaw', 'pAdj', 'Sig');
    fprintf('%s\n', repmat('-', 1, 72));
    for r = 1:nPairs
        fprintf('%-22s  %-22s  %8.4g  %8.4g  %s\n', ...
            comboA{r}, comboB{r}, rawP(r), adjP(r), sigLabel{r});
    end
    fprintf('\n  Significant pairs: %d / %d\n\n', numel(sigPairs), nPairs);
end

end   % ── end main function ──────────────────────────────────────────────


% =========================================================================
% LOCAL HELPERS
% =========================================================================

function pairs = allPairs(keys)
%ALLPAIRS  All unique unordered pairs from a cell array of strings.
    n     = numel(keys);
    pairs = {};
    for i = 1:n
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