function [intTbl, contrastTbl, lmeModel, workTbl] = runLMMInteractionStat(tbl, metricVec, options)
%RUNLMMINTERACTIONSTAT  Factorial LMM with genotype x condition interaction.
%
%   Fits:  metric ~ genotype * condition + (1|bioRep)
%
%   Use it to test directly whether the effect of a treatment (stiffness or
%   ROCK silencing) DIFFERS between genotypes (e.g. CT vs CCM2), instead of
%   comparing "significant in one group, not significant in the other".
%
%   Companion to runLMMStat. runLMMStat is NOT modified, so existing results
%   stay identical.
%
% -------------------------------------------------------------------------
% EXAMPLES (see also main_CombineAndPlot)
%
%   % Genotype x stiffness (gel experiments, Fig 4)
%   [intTbl, ctrTbl] = runLMMInteractionStat(T_Vin, T_Vin.FRETav, ...
%       'FactorLevels', {'1kPa','10kPa'});
%
%   % Genotype x ROCK (Fig 3). 'rest' is the reference level.
%   [intTbl, ctrTbl] = runLMMInteractionStat(T_Vin, T_Vin.FRETav, ...
%       'FactorLevels', {'rest','R1','R2'});
%
%   % VE-cadherin: same call with T_VECad and a VECad metric
%   [intTbl, ctrTbl] = runLMMInteractionStat(T_VECad, T_VECad.FRETav, ...
%       'FactorLevels', {'1kPa','10kPa'});
%
% -------------------------------------------------------------------------
% REQUIRED INPUTS
%   tbl        : table with the data (same table you give runLMMStat).
%   metricVec  : numeric vector, one value per row of tbl.
%
% OPTIONAL NAME-VALUE INPUTS
%   'GenotypeCol'    : column with genotype.            Default "diseaseCat"
%   'FactorCol'      : column with treatment/condition. Default "extrainfoCat"
%   'BioRepCol'      : biological replicate column.     Default "expCat"
%   'GenotypeLevels' : ordered levels, FIRST = reference. Default {'CT','CCM'}
%   'FactorLevels'   : ordered levels, FIRST = reference. Default: as found
%   'FilterCol'      : extra column to filter on.       Default "sampletypeCat"
%   'FilterLevel'    : level to keep in FilterCol.      Default "TS"
%                      (set FilterCol = "" to disable)
%   'CollapseRep'    : average technical reps per bio-rep (as in runLMMStat).
%                      Default true
%   'DFMethod'       : 'satterthwaite' (default) | 'residual'
%   'Alpha'          : significance threshold. Default 0.05
%   'Verbose'        : print results. Default true
%
% OUTPUTS
%   intTbl      : ANOVA table of the model terms. The row 'genotype:condition'
%                 is the interaction test: THE evidence for claims such as
%                 "more pronounced in CCM2" or "only in control cells".
%   contrastTbl : follow-up contrasts, FDR-corrected within each Type:
%                   'Simple effect'   : treatment effect within each genotype
%                   'Genotype effect' : CCM vs CT within each condition
%                   'Interaction'     : difference between the two treatment
%                                       effects (CCM effect minus CT effect)
%   lmeModel    : fitted LinearMixedModel.
%   workTbl     : table actually passed to fitlme.

% =========================================================================
% 1. Arguments
% =========================================================================
arguments
    tbl                         table
    metricVec                   (:,1) double

    options.GenotypeCol         string  = "diseaseCat"
    options.FactorCol           string  = "extrainfoCat"
    options.BioRepCol           string  = "expCat"
    options.GenotypeLevels      cell    = {'CT','CCM'}
    options.FactorLevels        cell    = {}
    options.FilterCol           string  = "sampletypeCat"
    options.FilterLevel         string  = "TS"
    options.CollapseRep         (1,1) logical = true
    options.DFMethod            string  = "satterthwaite"
    options.Alpha               (1,1) double ...
        {mustBeGreaterThan(options.Alpha,0), mustBeLessThan(options.Alpha,1)} = 0.05
    options.Verbose             (1,1) logical = true
end

gCol     = char(options.GenotypeCol);
fCol     = char(options.FactorCol);
bCol     = char(options.BioRepCol);
alphaVal = options.Alpha;
dfm      = char(options.DFMethod);

assert(height(tbl) == numel(metricVec), ...
    'runLMMInteractionStat: tbl height (%d) must match numel(metricVec) (%d).', ...
    height(tbl), numel(metricVec));
for c = {gCol, fCol, bCol}
    assert(ismember(c{1}, tbl.Properties.VariableNames), ...
        'runLMMInteractionStat: column "%s" not found in tbl.', c{1});
end

% =========================================================================
% 2. Select rows
% =========================================================================
gVec = string(tbl.(gCol));
fVec = string(tbl.(fCol));
bVec = string(tbl.(bCol));

keep = ~isnan(metricVec) & gVec ~= "" & fVec ~= "" & ~ismissing(bVec);

% optional filter (default: TS only, so tailless controls are excluded)
if strlength(options.FilterCol) > 0 && ...
        ismember(char(options.FilterCol), tbl.Properties.VariableNames)
    keep = keep & string(tbl.(char(options.FilterCol))) == options.FilterLevel;
end

gLevels = cellstr(options.GenotypeLevels);
if isempty(options.FactorLevels)
    fLevels = cellstr(unique(fVec(keep), 'stable'));
else
    fLevels = cellstr(options.FactorLevels);
end

keep = keep & ismember(gVec, gLevels) & ismember(fVec, fLevels);

nG = numel(gLevels);
nC = numel(fLevels);
assert(nG >= 2 && nC >= 2, ...
    'runLMMInteractionStat: need at least 2 genotype levels and 2 condition levels.');

% =========================================================================
% 3. Working table
% =========================================================================
workTbl           = table();
workTbl.metric    = metricVec(keep);
workTbl.genotype  = categorical(gVec(keep), gLevels, 'Ordinal', false);
workTbl.condition = categorical(fVec(keep), fLevels, 'Ordinal', false);
workTbl.bioRep    = categorical(bVec(keep));

if height(workTbl) == 0
    warning('runLMMInteractionStat: no valid rows after filtering.');
    [intTbl, contrastTbl, lmeModel] = deal(table(), table(), []);
    return
end

if options.CollapseRep
    % Average technical replicates per bio-rep x genotype x condition
    [gid, idG, idC, idB] = findgroups(workTbl.genotype, workTbl.condition, workTbl.bioRep);
    mVals = splitapply(@(v) mean(v, 'omitnan'), workTbl.metric, gid);
    workTbl           = table();
    workTbl.metric    = mVals;
    workTbl.genotype  = categorical(string(idG), gLevels, 'Ordinal', false);
    workTbl.condition = categorical(string(idC), fLevels, 'Ordinal', false);
    workTbl.bioRep    = idB;
    if options.Verbose
        fprintf('  After aggregation to bioRep means: %d rows\n', height(workTbl));
    end
end

% Which genotype x condition cells contain data?
hasData = false(nG, nC);
for ig = 1:nG
    for ic = 1:nC
        hasData(ig, ic) = any(workTbl.genotype == gLevels{ig} & ...
                              workTbl.condition == fLevels{ic});
    end
end
if ~all(hasData(:))
    [mg, mc] = find(~hasData);
    for k = 1:numel(mg)
        warning('runLMMInteractionStat: no data for %s x %s; the model is not fully crossed.', ...
            gLevels{mg(k)}, fLevels{mc(k)});
    end
end
assert(sum(hasData(:)) >= 4, ...
    'runLMMInteractionStat: too few filled cells to estimate an interaction.');

% =========================================================================
% 4. Fit the model
% =========================================================================
if options.Verbose
    fprintf('\n=== runLMMInteractionStat: fitting LMM ===\n');
    fprintf('  metric ~ genotype * condition + (1|bioRep)\n');
    fprintf('  Genotypes  (%d): %s  (reference: %s)\n', nG, strjoin(gLevels, ', '), gLevels{1});
    fprintf('  Conditions (%d): %s  (reference: %s)\n', nC, strjoin(fLevels, ', '), fLevels{1});
    fprintf('  Valid obs      : %d\n', height(workTbl));
    fprintf('  Bio replicates : %d\n', numel(unique(workTbl.bioRep)));
end

try
    lmeModel = fitlme(workTbl, 'metric ~ genotype*condition + (1|bioRep)', ...
        'FitMethod',      'REML', ...
        'DummyVarCoding', 'effects');
catch ME
    error('runLMMInteractionStat: fitlme failed: %s', ME.message);
end

% =========================================================================
% 5. Omnibus tests (main effects + interaction)
% =========================================================================
A = anova(lmeModel, 'DFMethod', dfm);
termNames = cellstr(A.Term);
intTbl = table(termNames, A.FStat, A.DF1, A.DF2, A.pValue, ...
    'VariableNames', {'Term','F','DF1','DF2','pValue'});
intTbl.Sig = arrayfun(@starsFromP, intTbl.pValue, 'UniformOutput', false);

% =========================================================================
% 6. Follow-up contrasts (cell-mean vectors from the design matrix)
% =========================================================================
X   = designMatrix(lmeModel, 'Fixed');
nFE = size(X, 2);

cm = zeros(nG, nC, nFE);
for ig = 1:nG
    for ic = 1:nC
        rows = workTbl.genotype == gLevels{ig} & workTbl.condition == fLevels{ic};
        if any(rows)
            cm(ig, ic, :) = mean(X(rows, :), 1);
        end
    end
end
getCM = @(ig, ic) reshape(cm(ig, ic, :), 1, nFE);

betaHat = fixedEffects(lmeModel);
CovBeta = lmeModel.CoefficientCovariance;

C = struct('Type',{}, 'Label',{}, 'L',{});

% (a) Simple effects: treatment effect within each genotype
for ig = 1:nG
    for i = 1:nC
        for j = i+1:nC
            if hasData(ig, i) && hasData(ig, j)
                C(end+1) = struct('Type', 'Simple effect', ...
                    'Label', sprintf('%s: %s vs %s', gLevels{ig}, fLevels{j}, fLevels{i}), ...
                    'L', getCM(ig, j) - getCM(ig, i)); %#ok<AGROW>
            end
        end
    end
end

% (b) Genotype effect within each condition
for ic = 1:nC
    for g1 = 1:nG
        for g2 = g1+1:nG
            if hasData(g1, ic) && hasData(g2, ic)
                C(end+1) = struct('Type', 'Genotype effect', ...
                    'Label', sprintf('%s: %s vs %s', fLevels{ic}, gLevels{g2}, gLevels{g1}), ...
                    'L', getCM(g2, ic) - getCM(g1, ic)); %#ok<AGROW>
            end
        end
    end
end

% (c) Interaction contrasts: (effect in genotype 2) - (effect in genotype 1)
for g1 = 1:nG
    for g2 = g1+1:nG
        for i = 1:nC
            for j = i+1:nC
                if all([hasData(g1,i), hasData(g1,j), hasData(g2,i), hasData(g2,j)])
                    C(end+1) = struct('Type', 'Interaction', ...
                        'Label', sprintf('(%s vs %s) in %s minus in %s', ...
                                fLevels{j}, fLevels{i}, gLevels{g2}, gLevels{g1}), ...
                        'L', (getCM(g2,j) - getCM(g2,i)) - (getCM(g1,j) - getCM(g1,i))); %#ok<AGROW>
                end
            end
        end
    end
end

nC_tot = numel(C);
est = nan(nC_tot,1); se = nan(nC_tot,1); tS = nan(nC_tot,1);
dfv = nan(nC_tot,1); pRaw = nan(nC_tot,1);

for r = 1:nC_tot
    L = C(r).L;
    [p, F, ~, df2] = coefTest(lmeModel, L, 0, 'DFMethod', dfm);
    est(r)  = L * betaHat;
    se(r)   = sqrt(max(0, L * CovBeta * L'));
    tS(r)   = est(r) / se(r);
    dfv(r)  = df2;
    pRaw(r) = p;
    %#ok<NASGU> F is unused (t = estimate/SE is reported instead)
end

% FDR (Benjamini-Hochberg) within each Type of contrast
typeVec = string({C.Type})';
pAdj = nan(nC_tot,1);
for t = unique(typeVec)'
    m = typeVec == t;
    pAdj(m) = bhAdjust(pRaw(m));
end

contrastTbl = table(typeVec, string({C.Label})', est, se, tS, dfv, pRaw, pAdj, ...
    arrayfun(@starsFromP, pAdj, 'UniformOutput', false), ...
    'VariableNames', {'Type','Contrast','Estimate','SE','tStat','DF','pRaw','pAdj','Sig'});

% =========================================================================
% 7. Verbose output
% =========================================================================
if options.Verbose
    disp(lmeModel);
    fprintf('\n--- Omnibus tests (DF method: %s) ---\n', dfm);
    fprintf('%-24s %8s %6s %8s %10s  %s\n', 'Term','F','DF1','DF2','p','Sig');
    for r = 1:height(intTbl)
        fprintf('%-24s %8.3f %6.0f %8.1f %10.4g  %s\n', intTbl.Term{r}, ...
            intTbl.F(r), intTbl.DF1(r), intTbl.DF2(r), intTbl.pValue(r), intTbl.Sig{r});
    end
    fprintf('  >> The ''genotype:condition'' row is the interaction test.\n');

    fprintf('\n--- Follow-up contrasts (FDR within each Type, alpha = %.2f) ---\n', alphaVal);
    for t = unique(typeVec, 'stable')'
        fprintf('\n[%s]\n', t);
        idx = find(typeVec == t);
        for r = idx'
            fprintf('  %-46s est=%9.4g  pRaw=%8.4g  pAdj=%8.4g  %s\n', ...
                C(r).Label, est(r), pRaw(r), pAdj(r), contrastTbl.Sig{r});
        end
    end
    fprintf('\n');
end
end


% =========================================================================
% LOCAL HELPERS
% =========================================================================
function adj = bhAdjust(p)
% Benjamini-Hochberg adjusted p-values (no toolbox needed)
    p = p(:);
    n = numel(p);
    [ps, ord] = sort(p);
    a = ps .* n ./ (1:n)';
    a = min(1, flipud(cummin(flipud(a))));
    adj = nan(n,1);
    adj(ord) = a;
end

function out = starsFromP(p)
    if isnan(p),      out = 'NA';
    elseif p < 0.001, out = '***';
    elseif p < 0.01,  out = '**';
    elseif p < 0.05,  out = '*';
    else,             out = 'ns';
    end
end