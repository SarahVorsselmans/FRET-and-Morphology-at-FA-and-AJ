function [sigPairs, statsTbl, mdl, workTbl] = runWeightedContrastTest(diffVec, sdVec, groupVec, options)
arguments
    diffVec   (:,1) double
    sdVec     (:,1) double
    groupVec
    options.Contrasts  = 'all'
    options.Correction string = "fdr"
    options.Alpha      (1,1) double = 0.05
    options.Verbose    (1,1) logical = true
end

groupVec = string(groupVec);
validMask = ~isnan(diffVec) & ~isnan(sdVec) & sdVec > 0;

workTbl        = table();
workTbl.diff   = diffVec(validMask);
workTbl.sd     = sdVec(validMask);
workTbl.group  = categorical(groupVec(validMask));
workTbl.weight = 1 ./ (workTbl.sd .^ 2);

groupLevels = categories(workTbl.group);
nGroup = numel(groupLevels);
if nGroup < 2
    warning('Fewer than 2 groups; nothing to test.');
    [sigPairs, statsTbl, mdl] = deal({}, table(), []);
    return
end

mdl = fitlm(workTbl, 'diff ~ group', 'Weights', workTbl.weight);
if options.Verbose, disp(mdl); end

if strcmpi(string(options.Contrasts), "all")
    pairs = {};
    for i = 1:nGroup
        for j = i+1:nGroup
            pairs{end+1} = {groupLevels{i}, groupLevels{j}}; %#ok<AGROW>
        end
    end
else
    pairs = options.Contrasts;
end

nPairs = numel(pairs);
rawP = nan(nPairs,1); est = nan(nPairs,1); se = nan(nPairs,1);

coefNames = mdl.CoefficientNames;   % e.g. {'(Intercept)', 'group_CCM_...'}
nCoef     = mdl.NumCoefficients;

% Which level is the reference? The one with no matching coefficient.
refGroup = groupLevels(~ismember( ...
    "group_" + string(groupLevels), string(coefNames)));
if ~isempty(refGroup), refGroup = refGroup{1}; else, refGroup = ""; end

for r = 1:nPairs

    gA = pairs{r}{1}; gB = pairs{r}{2};
    L = zeros(1, nCoef);

    if ~strcmp(gA, refGroup)
        idxA = strcmp(coefNames, "group_" + gA);
        L(idxA) = L(idxA) + 1;
    end
    if ~strcmp(gB, refGroup)
        idxB = strcmp(coefNames, "group_" + gB);
        L(idxB) = L(idxB) - 1;
    end

    [p, ~, ~] = coefTest(mdl, L);
    rawP(r) = p;
    est(r)  = L * mdl.Coefficients.Estimate;
    se(r)   = sqrt(L * mdl.CoefficientCovariance * L');
end

switch lower(char(options.Correction))
    case 'fdr',        adjP = mafdr(rawP, 'BHFDR', true);
    case 'bonferroni', adjP = min(rawP * nPairs, 1);
    otherwise,         adjP = rawP;
end

comboA = cellfun(@(p) p{1}, pairs, 'UniformOutput', false)';
comboB = cellfun(@(p) p{2}, pairs, 'UniformOutput', false)';
sigLabel = arrayfun(@pToStars, adjP, 'UniformOutput', false);

statsTbl = table(comboA, comboB, est, se, rawP, adjP, sigLabel, ...
    'VariableNames', {'ComboA','ComboB','Estimate','SE','pRaw','pAdj','Sig'});

sigPairs = {};
for r = 1:nPairs
    if adjP(r) < options.Alpha
        sigPairs{end+1} = {comboA{r}, comboB{r}, sigLabel{r}}; %#ok<AGROW>
    end
end
end

function out = pToStars(p)
    if p < 0.001,     out = '***';
    elseif p < 0.01,  out = '**';
    elseif p < 0.05,  out = '*';
    else,             out = 'ns';
    end
end