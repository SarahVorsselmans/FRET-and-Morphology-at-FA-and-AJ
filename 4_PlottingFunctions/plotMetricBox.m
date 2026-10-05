function [fig, pVals, overviewTbl] = plotMetricBox(tableT, variable, levels, compareCat, xCat, varargin)
%PLOTMETRICSCATTER  Boxplot of a metric per x-category, split by two
%                   comparison groups (e.g. CT vs CCM, or TS vs TL).
%                   compare1 boxes are blue, compare2 boxes are red.
%                   A gray vertical band highlights the compare2 side.
%
% USAGE
%   [fig, pVals, overviewTbl] = plotMetricScatter( ...
%       tableT, variable, levels, compareCat, xCat, Name, Value, ...)
%
% REQUIRED INPUTS
%   tableT      : MATLAB table.
%   variable    : char/string. Name of the numeric metric column to plot.
%   levels      : cell array of strings. Ordered list of x-axis categories.
%   compareCat  : char/string. Column name defining the two groups.
%   xCat        : char/string. Column name driving the x-axis grouping.
%
% OPTIONAL NAME-VALUE INPUTS
%   'compare1'        : Value in compareCat for the LEFT group (control).
%   'compare2'        : Value in compareCat for the RIGHT group (disease/test).
%   'includePatterns' : Cell array of strings. Default: {} (keep all).
%   'excludePatterns' : Cell array of strings. Default: {}.
%   'cutoff_bottom'   : Numeric scalar. Default: -Inf.
%   'cutoff_top'      : Numeric scalar. Default: Inf.
%   'splitCat'        : char/string. Column name used to pre-filter rows.
%   'splitCatLevel'   : string/char. Value to keep in splitCat column.
%   'StatTest'        : 'welch' (default) or 'mannwhitney'.
%   'BoxWidth'        : Scalar width of each box. Default: 0.12.
%
% OUTPUTS
%   fig         : Figure handle (or [] if no valid data).
%   pVals       : 1-by-K vector of p-values (one per level).
%   overviewTbl : Table with mean, SD, N per group per level, plus p-values.

% =========================================================================
% 1. Parse name-value inputs
% =========================================================================
p = inputParser;
p.KeepUnmatched = false;

addRequired(p, 'tableT',     @istable);
addRequired(p, 'variable',   @(x) ischar(x)||isstring(x));
addRequired(p, 'levels',     @(x) iscell(x)||isstring(x));
addRequired(p, 'compareCat', @(x) ischar(x)||isstring(x));
addRequired(p, 'xCat',       @(x) ischar(x)||isstring(x));

addParameter(p, 'compare1',        '',    @(x) ischar(x)||isstring(x)||isnumeric(x));
addParameter(p, 'compare2',        '',    @(x) ischar(x)||isstring(x)||isnumeric(x));
addParameter(p, 'includePatterns', {},    @iscell);
addParameter(p, 'excludePatterns', {},    @iscell);
addParameter(p, 'cutoff_bottom',   -Inf,  @(x) isnumeric(x)&&isscalar(x));
addParameter(p, 'cutoff_top',       Inf,  @(x) isnumeric(x)&&isscalar(x));
addParameter(p, 'splitCat',        '',    @(x) ischar(x)||isstring(x));
addParameter(p, 'splitCatLevel',   '',    @(x) ischar(x)||isstring(x));
addParameter(p, 'StatTest',        'welch', @(x) ischar(x)||isstring(x));
addParameter(p, 'BoxWidth',        0.12,  @(x) isnumeric(x)&&isscalar(x)&&x>0);

parse(p, tableT, variable, levels, compareCat, xCat, varargin{:});
r = p.Results;

% Normalise string inputs
varName     = char(string(r.variable));
compareCat  = char(string(r.compareCat));
xCat        = char(string(r.xCat));
levels      = cellstr(string(r.levels));
statTest    = lower(char(string(r.StatTest)));
boxWidth    = r.BoxWidth;

% Colors
color1 = [0.22 0.45 0.75];   % blue  – compare1
color2 = [0.85 0.20 0.20];   % red   – compare2

% =========================================================================
% 2. Filter rows via include/exclude on SampleFolder
% =========================================================================
allData = Filter_IncludeExclude(tableT, r.includePatterns, r.excludePatterns);

% =========================================================================
% 3. Optional pre-filter: keep only rows where splitCat == splitCatLevel
% =========================================================================
splitCat      = char(string(r.splitCat));
splitCatLevel = char(string(r.splitCatLevel));

if ~isempty(splitCat)
    if ~ismember(splitCat, allData.Properties.VariableNames)
        error('splitCat column "%s" not found in table.', splitCat);
    end
    if isempty(splitCatLevel)
        error('splitCatLevel must be specified when splitCat is used.');
    end
    colSplit = string(allData.(splitCat));
    allData  = allData(colSplit == string(splitCatLevel), :);
    if height(allData) == 0
        warning('No rows remain after splitCat filter "%s" == "%s".', splitCat, splitCatLevel);
        fig = []; pVals = []; overviewTbl = table();
        return;
    end
end

% =========================================================================
% 4. Validate required columns
% =========================================================================
requiredCols = {compareCat, xCat, 'DataFolder', varName};
missing = setdiff(requiredCols, allData.Properties.VariableNames);
if ~isempty(missing)
    error('Missing required column(s): %s', strjoin(missing, ', '));
end
if ~isnumeric(allData.(varName))
    error('Column "%s" must be numeric.', varName);
end

% =========================================================================
% 5. Determine compare1 / compare2 automatically if not supplied
% =========================================================================
compareVals = unique(string(allData.(compareCat)), 'stable');
compareVals = compareVals(compareVals ~= "");

compare1 = string(r.compare1);
compare2 = string(r.compare2);

if compare1 == ""
    if numel(compareVals) < 1
        error('No values found in compareCat column "%s".', compareCat);
    end
    compare1 = compareVals(1);
    fprintf('compare1 auto-set to "%s"\n', compare1);
end
if compare2 == ""
    if numel(compareVals) < 2
        error('Only one unique value in compareCat "%s"; cannot auto-set compare2.', compareCat);
    end
    compare2 = compareVals(2);
    fprintf('compare2 auto-set to "%s"\n', compare2);
end

% =========================================================================
% 6. Apply cutoffs and build filter mask
% =========================================================================
metric     = double(allData.(varName));
metric(metric <  r.cutoff_bottom) = NaN;
metric(metric >  r.cutoff_top)    = NaN;

compareStr = string(allData.(compareCat));
isGrp1     = (compareStr == compare1);
isGrp2     = (compareStr == compare2);
isEither   = isGrp1 | isGrp2;

xCatCol = string(allData.(xCat));
xCatCat = categorical(xCatCol, levels, 'Ordinal', true);

validRows = ~isnan(metric) & isEither & ~ismissing(xCatCat);

metricF  = metric(validRows);
compareF = compareStr(validRows);
xCatF    = xCatCat(validRows);

if isempty(metricF)
    warning('No valid "%s" entries after filtering.', varName);
    fig = []; pVals = []; overviewTbl = table();
    return;
end

% =========================================================================
% 7. Levels actually present in filtered data
% =========================================================================
levelsPresent = levels(ismember(levels, cellstr(unique(xCatF))));
K = numel(levelsPresent);
if K == 0
    warning('No xCat levels present after filtering.');
    fig = []; pVals = []; overviewTbl = table();
    return;
end

% =========================================================================
% 8. Layout constants
% =========================================================================
dx      = 0.18;   % half-separation between the two boxes at each tick
bandW   = 0.15;   % half-width of gray band around compare2 centre

% =========================================================================
% 9. Pre-compute stats
% =========================================================================
pVals     = nan(1, K);
starStrs  = strings(1, K);
yMaxLocal = nan(1, K);
nGrp1     = zeros(1, K);
nGrp2     = zeros(1, K);

for k = 1:K
    lvl  = levelsPresent{k};
    idx1 = (compareF == compare1) & (xCatF == lvl);
    idx2 = (compareF == compare2) & (xCatF == lvl);
    g1   = metricF(idx1);  g1 = g1(~isnan(g1));
    g2   = metricF(idx2);  g2 = g2(~isnan(g2));
    nGrp1(k) = numel(g1);
    nGrp2(k) = numel(g2);

    if numel(g1) >= 2 && numel(g2) >= 2
        [pVals(k), starStrs(k), ~] = runGroupStat(g1, g2, statTest);
        yMaxLocal(k) = max([g1; g2]);
    end
end

% =========================================================================
% 10. Draw figure
% =========================================================================
fig = figure('Color', 'w');
ax  = axes('Parent', fig);
hold(ax, 'on');

% --- Gray bands for compare2 (drawn first) ---
for k = 1:K
    xLeft  = (k + dx) - bandW;
    xRight = (k + dx) + bandW;
    patch(ax, [xLeft xRight xRight xLeft], [0 0 1 1], ...
        [0.88 0.88 0.88], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.55, ...
        'Tag', sprintf('band_%d', k));
end

% =========================================================================
% 11. Draw boxplots
% =========================================================================
% Dummy handles for legend (invisible scatter, just to capture colors)
hLeg1 = patch(ax, NaN, NaN, color1, 'EdgeColor', color1*0.6, 'DisplayName', char(compare1));
hLeg2 = patch(ax, NaN, NaN, color2, 'EdgeColor', color2*0.6, 'DisplayName', char(compare2));

for k = 1:K
    lvl = levelsPresent{k};

    % --- compare1 (left, blue) ---
    g1 = metricF((compareF == compare1) & (xCatF == lvl));
    g1 = g1(~isnan(g1));
    if numel(g1) >= 1
        drawBox(ax, g1, k - dx, boxWidth, color1);
    end

    % --- compare2 (right, red) ---
    g2 = metricF((compareF == compare2) & (xCatF == lvl));
    g2 = g2(~isnan(g2));
    if numel(g2) >= 1
        drawBox(ax, g2, k + dx, boxWidth, color2);
    end
end

% =========================================================================
% 12. Fix x-limits and stretch gray bands to full y-extent
% =========================================================================
xlim(ax, [0.5, K + 0.5]);
set(fig, 'Renderer', 'opengl');

validMetric = metricF(~isnan(metricF));
padAll = 0.08 * range(validMetric);
if ~isfinite(padAll) || padAll == 0
    padAll = 0.05 * max(abs(validMetric));
    if ~isfinite(padAll), padAll = 1; end
end
yl = ylim(ax);
ylim(ax, [yl(1), yl(2) + padAll]);
yl = ylim(ax);

% Stretch gray bands to full y-range
for k = 1:K
    b = findobj(ax, 'Tag', sprintf('band_%d', k));
    if isempty(b), continue; end
    set(b, 'YData', [yl(1) yl(1) yl(2) yl(2)]);
end

% Push bands behind boxes
bands = findobj(ax, '-regexp', 'Tag', '^band_');
uistack(bands, 'bottom');

% =========================================================================
% 13. Significance bars and stars
% =========================================================================
for k = 1:K
    if strlength(starStrs(k)) > 0 && isfinite(yMaxLocal(k))
        yl  = ylim(ax);
        pad = max(0.03 * (yl(2) - yl(1)), 1e-12);
        yPos = yMaxLocal(k) + pad;

        hl = plot(ax, [k-dx, k+dx], [yPos, yPos], 'k-', 'LineWidth', 1.3);
        ht = text(ax, k, yPos + 0.01*(yl(2)-yl(1)), starStrs(k), ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k');

        uistack(hl, 'top');
        uistack(ht, 'top');

        yl = ylim(ax);
        if yPos * 1.05 > yl(2)
            ylim(ax, [yl(1), yPos * 1.10]);
        end
    end
end

% =========================================================================
% 14. N labels above plot
% =========================================================================
yl = ylim(ax);
nLabelPad = 0.01 * (yl(2) - yl(1));

for k = 1:K
    yl   = ylim(ax);
    yTop = yl(2) + nLabelPad;

    if nGrp1(k) > 0
        text(ax, k - dx, yTop, sprintf('n=%d', nGrp1(k)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'bottom', ...
            'FontSize',            7, ...
            'Color',               color1, ...
            'Clipping',            'off');
    end
    if nGrp2(k) > 0
        text(ax, k + dx, yTop, sprintf('n=%d', nGrp2(k)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'bottom', ...
            'FontSize',            7, ...
            'Color',               color2, ...
            'Clipping',            'off');
    end
end

yl = ylim(ax);
ylim(ax, [yl(1), yl(2) + 0.05 * (yl(2) - yl(1))]);

% =========================================================================
% 15. Axes decoration
% =========================================================================
xticks(ax, 1:K);
xticklabels(ax, levelsPresent);
ylabel(ax, varName, 'Interpreter', 'none');

titleStr = sprintf('%s  |  %s vs %s', varName, compare1, compare2);
title(ax, titleStr, 'Interpreter', 'none');

subtitleParts = {sprintf('stat: %s | gray band = %s', statTest, compare2)};
if isfinite(r.cutoff_bottom), subtitleParts{end+1} = sprintf('cut_bot=%.4g', r.cutoff_bottom); end
if isfinite(r.cutoff_top),    subtitleParts{end+1} = sprintf('cut_top=%.4g', r.cutoff_top);    end
subtitle(ax, strjoin(subtitleParts, '  |  '), 'FontSize', 8, 'Color', [0.4 0.4 0.4]);

grid(ax, 'on');
set(ax, 'Box', 'off', 'Layer', 'top');

% Legend: blue = compare1, red = compare2
legend([hLeg1, hLeg2], ...
    {char(compare1), char(compare2)}, ...
    'Location',    'southoutside', ...
    'Interpreter', 'none', ...
    'NumColumns',  2);

% =========================================================================
% 16. Overview table
% =========================================================================
varTypes = {'string','double','double','double','double','double','double','double','double','string'};
varNames = {'Level', ...
            sprintf('Mean_%s', compare1), sprintf('SD_%s', compare1), sprintf('N_%s', compare1), ...
            sprintf('Mean_%s', compare2), sprintf('SD_%s', compare2), sprintf('N_%s', compare2), ...
            'Diff_2_minus_1', 'pValue', 'Significance'};

overviewTbl = table('Size', [K numel(varNames)], ...
    'VariableTypes', varTypes, ...
    'VariableNames', varNames);

for k = 1:K
    lvl  = levelsPresent{k};
    idx1 = (compareF == compare1) & (xCatF == lvl);
    idx2 = (compareF == compare2) & (xCatF == lvl);
    g1   = metricF(idx1);  g1 = g1(~isnan(g1));
    g2   = metricF(idx2);  g2 = g2(~isnan(g2));

    overviewTbl.Level(k)                            = string(lvl);
    overviewTbl.(sprintf('Mean_%s', compare1))(k)  = mean(g1, 'omitnan');
    overviewTbl.(sprintf('SD_%s',   compare1))(k)  = std(g1,  'omitnan');
    overviewTbl.(sprintf('N_%s',    compare1))(k)  = numel(g1);
    overviewTbl.(sprintf('Mean_%s', compare2))(k)  = mean(g2, 'omitnan');
    overviewTbl.(sprintf('SD_%s',   compare2))(k)  = std(g2,  'omitnan');
    overviewTbl.(sprintf('N_%s',    compare2))(k)  = numel(g2);
    overviewTbl.Diff_2_minus_1(k)                  = mean(g2,'omitnan') - mean(g1,'omitnan');
    overviewTbl.pValue(k)                          = pVals(k);
    overviewTbl.Significance(k)                    = string(starStrs(k));
end

fprintf('\n--- Overview: %s vs %s per %s ---\n', compare1, compare2, xCat);
disp(overviewTbl);

hold(ax, 'off');
end


% =========================================================================
% LOCAL HELPER: drawBox
%   Draws a single styled boxplot at position xCenter.
%
%   Inputs
%     ax       : target axes handle
%     data     : numeric vector (NaN-free)
%     xCenter  : scalar x position
%     bw       : full box width in x-axis units
%     clr      : 1x3 RGB face color
% =========================================================================
function drawBox(ax, data, xCenter, bw, clr)

    if isempty(data)
        return;
    end

    % Compute box statistics
    q1  = quantile(data, 0.25);
    med = median(data);
    q3  = quantile(data, 0.75);
    iqr = q3 - q1;

    % Whisker ends (Tukey style, capped at data extremes)
    wLow  = max(min(data), q1 - 1.5*iqr);
    wHigh = min(max(data), q3 + 1.5*iqr);

    % Outliers
    outliers = data(data < wLow | data > wHigh);

    hw = bw / 2;   % half box width

    % Box fill
    patch(ax, ...
        [xCenter-hw, xCenter+hw, xCenter+hw, xCenter-hw], ...
        [q1, q1, q3, q3], ...
        clr, ...
        'FaceAlpha', 0.45, ...
        'EdgeColor', clr * 0.6, ...
        'LineWidth', 1.2);

    % Median line
    plot(ax, [xCenter-hw, xCenter+hw], [med, med], ...
        'Color', clr * 0.5, 'LineWidth', 2.0);

    % Whiskers
    plot(ax, [xCenter, xCenter], [wLow,  q1], '-', 'Color', clr*0.6, 'LineWidth', 1.0);
    plot(ax, [xCenter, xCenter], [q3, wHigh], '-', 'Color', clr*0.6, 'LineWidth', 1.0);

    % Whisker caps
    capW = hw * 0.5;
    plot(ax, [xCenter-capW, xCenter+capW], [wLow,  wLow],  '-', 'Color', clr*0.6, 'LineWidth', 1.0);
    plot(ax, [xCenter-capW, xCenter+capW], [wHigh, wHigh], '-', 'Color', clr*0.6, 'LineWidth', 1.0);

    % Outliers
    if ~isempty(outliers)
        scatter(ax, repmat(xCenter, numel(outliers), 1), outliers, 18, ...
            'MarkerEdgeColor', clr*0.6, ...
            'MarkerFaceColor', clr, ...
            'MarkerFaceAlpha', 0.5);
    end
end