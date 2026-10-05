function [pTbl, anovaTbl, statsStruct] = runANOVAConditions(data, conditionLabels, groupLabels, varargin)
%RUNANOVACONDITIONS  Two-way ANOVA with post-hoc tests across conditions,
%                   split by two groups (e.g. control vs disease).
%
%   Runs a two-way ANOVA with factors Condition (xCat levels) and Group
%   (compare1 vs compare2), including their interaction. Then runs a
%   post-hoc multiple-comparison test across conditions — separately for
%   each group and combined — to identify which conditions differ.
%
% USAGE
%   [pTbl, anovaTbl, statsStruct] = runANOVAConditions( ...
%       data, conditionLabels, groupLabels, Name, Value)
%
% REQUIRED INPUTS
%   data            : numeric vector. All observations pooled.
%   conditionLabels : string/cell/categorical vector. Condition per obs.
%   groupLabels     : string/cell/categorical vector. Group per obs.
%
% OPTIONAL NAME-VALUE INPUTS
%   'PostHocTest'   : 'games-howell' (default) | 'tukey-kramer' | 'tukey'
%   'Alpha'         : Significance level. Default: 0.05.
%   'Verbose'       : Print results to Command Window. Default: true.
%
% OUTPUTS
%   pTbl        : Table with columns:
%                   Group, Condition_A, Condition_B, pValue, Significance
%   anovaTbl    : ANOVA table from anovan (cell array).
%   statsStruct : Struct with raw multcompare outputs per subset.

% =========================================================================
% 1. Parse inputs manually (no arguments block — handles categorical inputs)
% =========================================================================
p = inputParser;
p.KeepUnmatched = false;
addParameter(p, 'PostHocTest', 'games-howell', @(x) ischar(x) || isstring(x));
addParameter(p, 'Alpha',       0.05,           @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Verbose',     true,           @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});

postHocTest = lower(char(string(p.Results.PostHocTest)));
alpha       = p.Results.Alpha;
verbose     = logical(p.Results.Verbose);

% =========================================================================
% 2. Sanitise inputs — force correct shapes and types
% =========================================================================
data            = double(data(:));
conditionLabels = string(conditionLabels(:));
groupLabels     = string(groupLabels(:));

if numel(data) ~= numel(conditionLabels) || numel(data) ~= numel(groupLabels)
    error('runANOVAConditions: data, conditionLabels, and groupLabels must have the same length.');
end

% Remove NaNs
valid           = ~isnan(data);
data            = data(valid);
conditionLabels = conditionLabels(valid);
groupLabels     = groupLabels(valid);

if numel(data) < 4
    error('runANOVAConditions: fewer than 4 valid observations after NaN removal.');
end

% Unique levels
conditions = unique(conditionLabels, 'stable');
groups     = unique(groupLabels,     'stable');
nCond      = numel(conditions);

if nCond < 2
    error('runANOVAConditions: need at least 2 condition levels.');
end

% =========================================================================
% 3. Two-way ANOVA  (condition + group + interaction)
% =========================================================================
if verbose
    fprintf('\n=== Two-Way ANOVA: Condition x Group ===\n');
    fprintf('    Post-hoc: %s | alpha = %.3f\n\n', postHocTest, alpha);
end

% Convert to cell arrays of chars — most compatible anovan input
condCell  = cellstr(conditionLabels);
groupCell = cellstr(groupLabels);

[pAnova, anovaTbl, ~] = anovan(data, ...
    {condCell, groupCell}, ...
    'model',    'interaction', ...
    'varnames', {'Condition', 'Group'}, ...
    'display',  'off');

if verbose
    fprintf('ANOVA p-values:\n');
    fprintf('  Condition  (main effect)        : p = %.4g\n', pAnova(1));
    fprintf('  Group      (main effect)        : p = %.4g\n', pAnova(2));
    fprintf('  Condition x Group (interaction) : p = %.4g\n\n', pAnova(3));
end

% =========================================================================
% 4. Map post-hoc choice to multcompare CType
% =========================================================================
switch postHocTest
    case {'games-howell', 'tukey-kramer'}
        cType = 'hsd';           % Tukey-Kramer / HSD — handles unequal n
    case 'tukey'
        cType = 'tukey-kramer';  % MATLAB's name for balanced Tukey
    otherwise
        warning('runANOVAConditions: unknown PostHocTest "%s", using hsd.', postHocTest);
        cType = 'hsd';
end

% =========================================================================
% 5. Post-hoc pairwise condition comparisons
%    Run separately for: group1 only, group2 only
%    Plus within-level group1 vs group2 per condition pair
% =========================================================================
statsStruct  = struct();
pairRows     = {};

subsetNames  = cellstr(groups(:)');
subsetLabels = cellstr(groups(:)');

% --- Per-group across-level comparisons (unchanged logic) ---
for s = 1:numel(subsetNames)
    sName  = subsetNames{s};
    sLabel = subsetLabels{s};

    mask    = strcmp(groupLabels, sName);
    subData = data(mask);
    subCond = cellstr(conditionLabels(mask));

    if numel(unique(subCond)) < 2 || numel(subData) < 4
        if verbose
            fprintf('Skipping subset "%s": insufficient data.\n', sLabel);
        end
        continue;
    end

    [uCond, ~, idx] = unique(subCond);
    counts = accumarray(idx, 1);
    validCondMask = counts >= 2;
    if sum(validCondMask) < 2
        if verbose
            fprintf('Skipping subset "%s": not enough valid conditions.\n', sLabel);
        end
        continue;
    end
    validLevels = uCond(validCondMask);
    keepMask    = ismember(subCond, validLevels);
    subData     = subData(keepMask);
    subCond     = subCond(keepMask);

    [~, ~, statsOne] = anova1(subData, subCond, 'off');

    assert(ischar(cType) && isrow(cType), 'cType must be a row char vector')
    [compMatrix, means] = multcompare(statsOne, ...
        'CType',   cType, ...
        'Alpha',   alpha, ...
        'Display', 'off');

    safeField = matlab.lang.makeValidName(sName);
    statsStruct.(safeField).comparison = compMatrix;
    statsStruct.(safeField).means      = means;
    statsStruct.(safeField).stats      = statsOne;

    condList = statsOne.gnames;

    if verbose
        fprintf('--- Post-hoc (%s): %s across conditions ---\n', postHocTest, sLabel);
        fprintf('  %-20s  %-20s  %10s  %s\n', 'Condition A', 'Condition B', 'p-value', 'Sig');
    end

    for r = 1:size(compMatrix, 1)
        i     = compMatrix(r, 1);
        j     = compMatrix(r, 2);
        pVal  = compMatrix(r, 6);
        condA = condList{i};
        condB = condList{j};
        sig   = pToStars(pVal);

        if verbose
            fprintf('  %-20s  %-20s  %10.4g  %s\n', condA, condB, pVal, sig);
        end

        pairRows{end+1} = {string(sLabel), string(condA), string(condB), pVal, string(sig)};
    end

    if verbose, fprintf('\n'); end
end

% --- Within-level compare1 vs compare2 ---
if verbose
    fprintf('--- Within-level: %s vs %s per condition ---\n', groups(1), groups(2));
    fprintf('  %-20s  %10s  %s\n', 'Condition', 'p-value', 'Sig');
end

for c = 1:numel(conditions)
    lvl   = conditions(c);
    mask1 = strcmp(groupLabels, groups(1)) & strcmp(conditionLabels, lvl);
    mask2 = strcmp(groupLabels, groups(2)) & strcmp(conditionLabels, lvl);
    g1    = data(mask1);
    g2    = data(mask2);

    if numel(g1) < 2 || numel(g2) < 2
        if verbose
            fprintf('  %-20s  skipped (n<%d)\n', lvl, 2);
        end
        continue;
    end

    [~, pVal] = ttest2(g1, g2, 'Vartype', 'unequal');   % Welch t-test
    sig       = pToStars(pVal);

    if verbose
        fprintf('  %-20s  %10.4g  %s\n', lvl, pVal, sig);
    end

    pairRows{end+1} = {string(sprintf('%s vs %s', groups(1), groups(2))), ...
                       string(lvl), string(lvl), pVal, string(sig)};
end

if verbose, fprintf('\n'); end

% % =========================================================================
% % 5. Post-hoc pairwise condition comparisons
% %    Run separately for: combined, group1 only, group2 only
% % =========================================================================
% statsStruct  = struct();
% pairRows     = {};
% 
% subsetNames  = [{"combined"}, cellstr(groups(:)')];
% subsetLabels = [{"Combined"}, cellstr(groups(:)')];
% 
% for s = 1:numel(subsetNames)
%     sName  = subsetNames{s};
%     sLabel = subsetLabels{s};
% 
%     if strcmp(sName, 'combined')
%         mask = true(numel(data), 1);
%     else
%         mask = strcmp(groupLabels, sName);
%     end
% 
%     subData = data(mask);
%     subCond = cellstr(conditionLabels(mask));   % cell array of chars for anovan
% 
%     if numel(unique(subCond)) < 2 || numel(subData) < 4
%         if verbose
%             fprintf('Skipping subset "%s": insufficient data.\n', sLabel);
%         end
%         continue;
%     end
% 
%     % One-way ANOVA on this subset
%     [~, ~, statsOne] = anovan(subData, {subCond}, ...
%         'model',    'linear', ...
%         'varnames', {'Condition'}, ...
%         'display',  'off');
% 
%     % Post-hoc
%     [compMatrix, means] = multcompare(statsOne, ...
%         'CType',   cType, ...
%         'Alpha',   alpha, ...
%         'Display', 'off');
% 
%     % Store
%     safeField = matlab.lang.makeValidName(sName);
%     statsStruct.(safeField).comparison = compMatrix;
%     statsStruct.(safeField).means      = means;
%     statsStruct.(safeField).stats      = statsOne;
% 
%     % condList in the order anovan used (alphabetical by default)
%     condList = unique(subCond, 'sorted');
% 
%     if verbose
%         fprintf('--- Post-hoc (%s): %s across conditions ---\n', postHocTest, sLabel);
%         fprintf('  %-20s  %-20s  %10s  %s\n', 'Condition A', 'Condition B', 'p-value', 'Sig');
%     end
% 
%     for r = 1:size(compMatrix, 1)
%         i    = compMatrix(r, 1);
%         j    = compMatrix(r, 2);
%         pVal = compMatrix(r, 6);
%         condA = condList{i};
%         condB = condList{j};
%         sig   = pToStars(pVal);
% 
%         if verbose
%             fprintf('  %-20s  %-20s  %10.4g  %s\n', condA, condB, pVal, sig);
%         end
% 
%         pairRows{end+1} = {string(sLabel), string(condA), string(condB), pVal, string(sig)}; %#ok<AGROW>
%     end
% 
%     if verbose, fprintf('\n'); end
% end

% =========================================================================
% 6. Build output table
% =========================================================================
if isempty(pairRows)
    pTbl = table('Size', [0 5], ...
        'VariableTypes', {'string','string','string','double','string'}, ...
        'VariableNames', {'Group','Condition_A','Condition_B','pValue','Significance'});
else
    nRows = numel(pairRows);
    Group        = strings(nRows, 1);
    Condition_A  = strings(nRows, 1);
    Condition_B  = strings(nRows, 1);
    pValue       = nan(nRows, 1);
    Significance = strings(nRows, 1);

    for r = 1:nRows
        Group(r)        = pairRows{r}{1};
        Condition_A(r)  = pairRows{r}{2};
        Condition_B(r)  = pairRows{r}{3};
        pValue(r)       = pairRows{r}{4};
        Significance(r) = pairRows{r}{5};
    end

    pTbl = table(Group, Condition_A, Condition_B, pValue, Significance);
end

if verbose
    fprintf('=== Full pairwise table ===\n');
    disp(pTbl);
end

end % main function


% =========================================================================
% LOCAL HELPER: pToStars
% =========================================================================
function s = pToStars(p)
    if     p < 0.001,  s = '***';
    elseif p < 0.01,   s = '**';
    elseif p < 0.05,   s = '*';
    else,              s = 'ns';
    end
end
