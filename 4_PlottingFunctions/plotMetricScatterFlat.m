function [fig, anovaTbl, overviewTbl] = plotMetricScatterFlat(tableT, varName, options)
%PLOTMETRICSCATTERFLAT  Scatter plot of a metric with all (diseaseCat x
%                        extrainfoCat) group combinations laid out side by
%                        side on the x-axis.  No left/right paired split.
%
%                        Points are colored by DataFolder.  Group means are
%                        drawn as short horizontal lines.  Violin plots are
%                        drawn behind the scatter.  A one-way ANOVA with
%                        post-hoc testing is run across all groups, and
%                        significant pairwise brackets are drawn compactly
%                        above the data.
%
% USAGE
%   [fig, anovaTbl, overviewTbl] = plotMetricScatterFlat( ...
%       tableT, varName, Name, Value, ...)
%
% REQUIRED INPUTS
%   tableT    : MATLAB table.
%   varName   : char/string. Name of the numeric metric column to plot.
%
% OPTIONAL NAME-VALUE INPUTS
%   'GroupOrder'      : Cell array of strings {'disease_extra', ...}.
%                       Each entry is "<diseaseCat>_<extrainfoCat>" and
%                       defines x-position order.  If empty, all unique
%                       combinations found in the data are used in the
%                       order they appear.
%   'DiseaseCat'      : Column name for the disease/condition category.
%                       Default: "diseaseCat".
%   'XCat'            : Column name for the x-axis sub-category.
%                       Default: "extrainfoCat".
%   'includePatterns' : Cell array of strings. Keep rows matching any pattern.
%   'excludePatterns' : Cell array of strings. Drop rows matching any pattern.
%   'CutoffBottom'    : Numeric scalar. Values below are set to NaN.
%   'CutoffTop'       : Numeric scalar. Values above are set to NaN.
%   'YLimMaxSet'      : Override upper y-limit.
%   'YLimMinSet'      : Override lower y-limit.
%   'SplitCat'        : Column name for pre-filter.
%   'SplitCatLevel'   : Value to keep in SplitCat.
%   'PostHocTest'     : 'games-howell' (default) | 'tukey-kramer'.
%   'ANOVAVerbose'    : Logical. Print ANOVA table. Default: true.
%   'Palette'         : Color palette name. Default: 'distinct'.
%   'MarkerSize'      : Scalar. Default: 36.
%   'Alpha'           : Marker transparency [0,1]. Default: 0.75.
%   'Jitter'          : Horizontal jitter amplitude. Default: 0.08.
%   'ViolinWidth'     : Max half-width of each violin. Default: 0.13.
%   'ViolinAlpha'     : Violin transparency [0,1]. Default: 0.25.
%   'ShowViolin'      : Logical. Default: true.
%
% OUTPUTS
%   fig         : Figure handle (or [] if no valid data).
%   anovaTbl    : ANOVA table ([] if fewer than 2 groups present).
%   overviewTbl : Table with mean, SD, N per group.

% =========================================================================
% 1. Arguments
% =========================================================================
arguments
    tableT                table
    varName               {mustBeA(varName, ["string","char","double"])}
    options.VarColName    string = ""
    options.VarColIndex   (1,1) double {mustBePositive} = 1

    options.GroupOrder        cell   = {}
    options.DiseaseCat        string = "diseaseCat"
    options.XCat              string = "extrainfoCat"

    options.includePatterns   cell   = {}
    options.excludePatterns   cell   = {}

    options.CutoffBottom  (1,1) double = -Inf
    options.CutoffTop     (1,1) double = Inf
    options.YLimMaxSet    = []
    options.YLimMinSet    = []
    options.SplitCat      string = ""
    options.SplitCatLevel string = ""

    options.PostHocTest   string = "tukey-kramer"
    options.ANOVAVerbose  (1,1) logical = true

    options.paletteName   string = "distinct"
    options.MarkerSize    (1,1) double {mustBePositive} = 36
    options.Alpha         (1,1) double {mustBeGreaterThanOrEqual(options.Alpha,0),   mustBeLessThanOrEqual(options.Alpha,1)}   = 0.75
    options.Jitter        (1,1) double {mustBeNonnegative} = 0.08
    options.ViolinWidth   (1,1) double {mustBePositive} = 0.13
    options.ViolinAlpha   (1,1) double {mustBeGreaterThanOrEqual(options.ViolinAlpha,0), mustBeLessThanOrEqual(options.ViolinAlpha,1)} = 0.25
    options.ShowViolin    (1,1) logical = true

    options.VarNaming     string
    options.Unit          string = ""

    options.MultiplePalette      (1,1) logical = false
    options.MultiplePaletteNames cell  = {}
end

diseaseCat    = options.DiseaseCat;
xCat          = options.XCat;
groupOrder    = options.GroupOrder;
CutoffBottom  = options.CutoffBottom;
CutoffTop     = options.CutoffTop;
YLimMaxSet    = options.YLimMaxSet;
YLimMinSet    = options.YLimMinSet;
splitCat      = char(string(options.SplitCat));
splitCatLevel = char(string(options.SplitCatLevel));
postHocTest   = options.PostHocTest;
anovaVerbose  = options.ANOVAVerbose;
paletteName   = options.paletteName;
mSize         = options.MarkerSize;
alphaVal      = options.Alpha;
jitterAmp     = options.Jitter;
violinWidth   = options.ViolinWidth;
violinAlpha   = options.ViolinAlpha;
showViolin    = options.ShowViolin;
unitName      = options.Unit;

multiplePalette      = options.MultiplePalette;
multiplePaletteNames = options.MultiplePaletteNames;

% % if varName is not chr
% if ischar(varName) || isstring(varName)
%     metricRaw = tableT.(varName);
%     varNameStr = string(varName);
% else
%     metricRaw = varName;
%     varNameStr = options.VarNaming;
% end
% if varName is not chr
if ischar(varName) || isstring(varName)
    varNameStr = string(varName);
    metricIsName = true;
else
    % numeric vector case — requires VarColName + VarColIndex to re-extract
    varNameStr   = options.VarNaming;
    metricIsName = false;
    if options.VarColName == ""
        error('When passing a numeric vector as varName, you must also supply ''VarColName'' and ''VarColIndex''.');
    end
end

% =========================================================================
% 2. Include / exclude filter
% =========================================================================
allData = Filter_IncludeExclude(tableT, options.includePatterns, options.excludePatterns);

% =========================================================================
% 3. Optional SplitCat pre-filter
% =========================================================================
if ~isempty(splitCat)
    if ~ismember(splitCat, allData.Properties.VariableNames)
        error('SplitCat column "%s" not found in table.', splitCat);
    end
    if isempty(splitCatLevel)
        error('SplitCatLevel must be specified when SplitCat is used.');
    end
    allData = allData(string(allData.(splitCat)) == string(splitCatLevel), :);
    if height(allData) == 0
        warning('No rows remain after SplitCat filter "%s" == "%s".', splitCat, splitCatLevel);
        fig = []; anovaTbl = []; overviewTbl = table();
        return;
    end
end

% Extract metric from the already-filtered allData
if metricIsName
    metricRaw = allData.(varNameStr);
else
    metricRaw = allData.(options.VarColName)(:, options.VarColIndex);
end

% =========================================================================
% 4. Build combined group label per row: "<diseaseCat>_<extrainfoCat>"
% =========================================================================
if ~ismember(diseaseCat, allData.Properties.VariableNames)
    error('DiseaseCat column "%s" not found in table.', diseaseCat);
end
if ~ismember(xCat, allData.Properties.VariableNames)
    error('XCat column "%s" not found in table.', xCat);
end

combinedLabel = string(allData.(diseaseCat)) + "_" + string(allData.(xCat));

% =========================================================================
% 5. Determine group order
% =========================================================================
if isempty(groupOrder)
    % Preserve appearance order
    [~, ia]    = unique(combinedLabel, 'stable');
    groupOrder = cellstr(combinedLabel(sort(ia)));
end
G = numel(groupOrder);

% =========================================================================
% 6. Cutoffs and valid-row mask
% =========================================================================
%metricRaw = allData.(varName);
if iscell(metricRaw)
    metric = nan(size(metricRaw));
    for i = 1:numel(metricRaw)
        v = metricRaw{i};
        if ~isempty(v) && isnumeric(v)
            metric(i) = v;
        end
    end
else
    metric = double(metricRaw);
end
metric(metric < CutoffBottom) = NaN;
metric(metric > CutoffTop)    = NaN;

folder = categorical(string(allData.DataFolder));

groupCat  = categorical(combinedLabel, groupOrder, 'Ordinal', true);
validRows = ~isnan(metric) & ~ismissing(groupCat) & ~ismissing(folder);

metricF = metric(validRows);
groupF  = groupCat(validRows);
folderF = folder(validRows);

if isempty(metricF)
    warning('No valid "%s" entries after filtering.', varNameStr);
    fig = []; anovaTbl = []; overviewTbl = table();
    return;
end

% =========================================================================
% 7. Groups actually present after filtering
% =========================================================================
levelsPresent = groupOrder(ismember(groupOrder, cellstr(unique(groupF))));
K = numel(levelsPresent);
if K == 0
    warning('No groups present after filtering.');
    fig = []; anovaTbl = []; overviewTbl = table();
    return;
end

% =========================================================================
% 8. Colors and markers (by DataFolder)
% =========================================================================
folderCats = unique(folderF, 'stable');
nF         = numel(folderCats);
colorsF    = getPlotColors(paletteName, nF);
markersF   = getPlotMarkers(nF);

% =========================================================================
% 9. One-way ANOVA + post-hoc
% =========================================================================
% Issue with DIPImage interference -> 
dipPath = fileparts(which('split'));   % currently the DIPimage one
rmpath(dipPath); 

anovaTbl  = [];
sigPairs  = {};   % cell of {groupA, groupB, starStr} for significant pairs

if K >= 2
    groupVec = cellstr(string(groupF));

    % --- One-way ANOVA ---
    [pAnova, anovaTbl] = anova1(metricF, groupVec, 'off');

    if anovaVerbose
        fprintf('\n--- One-way ANOVA: %s ---\n', varNameStr);
        fprintf('  Overall p = %.4g\n', pAnova);
        disp(anovaTbl);
    end

    % --- Post-hoc pairwise comparison ---
    % multcompare needs the output of anova1 run with stats struct
    [~, ~, stats] = anova1(metricF, groupVec, 'off');

    try
        if strcmpi(postHocTest, 'games-howell')
            % Games-Howell is not in base MATLAB; fall back to Tukey-Kramer
            % if the Statistics Toolbox version doesn't support it.
            try
                [c, ~, ~, gnames] = multcompare(stats, ...
                    'CType', 'games-howell', 'Display', 'off');
            catch
                warning('games-howell not available; falling back to tukey-kramer.');
                [c, ~, ~, gnames] = multcompare(stats, ...
                    'CType', 'tukey-kramer', 'Display', 'off');
            end
        else
            [c, ~, ~, gnames] = multcompare(stats, ...
                'CType', postHocTest, 'Display', 'off');
        end

        % c columns: [i j lower diff upper p]
        for r = 1:size(c, 1)
            pPair = c(r, 6);
            if pPair < 0.05
                gA = char(gnames{c(r,1)});
                gB = char(gnames{c(r,2)});
                sigPairs{end+1} = {gA, gB, pToStars(pPair)}; %#ok<AGROW>
            end
        end

        if anovaVerbose
            fprintf('\n  Post-hoc (%s) significant pairs (p<0.05):\n', postHocTest);
            if isempty(sigPairs)
                fprintf('    none\n');
            else
                for r = 1:numel(sigPairs)
                    fprintf('    %s  vs  %s  %s\n', sigPairs{r}{1}, sigPairs{r}{2}, sigPairs{r}{3});
                end
            end
        end

    catch ME
        warning('%s', sprintf('Post-hoc test failed: %s', ME.message));
    end
end

addpath(dipPath); % restore DIPImage location

% =========================================================================
% 10. Layout constants
% =========================================================================
meanHalf = max(jitterAmp, 0.08);
dx       = 0.18;
xSpacing = 0.5;   % < 1 = tighter spacing
xPos = 1 + (0:K-1) * xSpacing;

% =========================================================================
% 11. Draw figure
% =========================================================================
fig = figure('Color', 'w');
ax  = axes('Parent', fig);
hold(ax, 'on');

% % --- Alternating light background bands for readability ---
% xPos = 1 + (0:K-1) * xSpacing;
%     if mod(k, 2) == 0
%         patch(ax, [k-0.5 k+0.5 k+0.5 k-0.5], [0 0 1 1], ...
%             [0.93 0.93 0.93], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
%             'Tag', sprintf('band_%d', k));
%     end
% end

% --- Violin plots ---
if showViolin
    violinColor = [0.60 0.60 0.70];
    nKde = 256;
    for k = 1:K
        grp = levelsPresent{k};
        g   = metricF(groupF == grp);
        g   = g(~isnan(g));
        if numel(g) >= 3
            drawViolin(ax, g, xPos(k), violinWidth, nKde, ...
                violinColor, violinAlpha, sprintf('violin_%d', k));
        end
    end
end


% --- MultiplePalette: per-group palettes ---
groupColorMap = containers.Map('KeyType','char','ValueType','any');
if multiplePalette
    if isempty(multiplePaletteNames)
        error('MultiplePaletteNames must be provided when MultiplePalette is true.');
    end
    for k = 1:K
        grp = levelsPresent{k};
        if k <= numel(multiplePaletteNames)
            palName = multiplePaletteNames{k};
        else
            warning('Not enough palette names supplied; recycling last one.');
            palName = multiplePaletteNames{end};
        end
        groupColorMap(char(grp)) = getPlotColors(palName, nF);
    end
end

% --- Scatter + mean lines ---
lgHandles = gobjects(nF, 1);
hasHandle = false(nF, 1);
nPerGroup = zeros(1, K);

for k = 1:K
    grp  = levelsPresent{k};
    idxG = (groupF == grp);
    nPerGroup(k) = sum(idxG & ~isnan(metricF));

    if ~any(idxG), continue; end

    yAll = metricF(idxG);
    for f = 1:nF
        idxF = idxG & (folderF == folderCats(f));
        yy   = metricF(idxF);
        if isempty(yy), continue; end
        xx = xPos(k) + jitterAmp * (rand(numel(yy), 1) - 0.5);
        if multiplePalette
            comboColors = groupColorMap(char(grp));
            dotColor = comboColors(f,:);
        else
            dotColor = colorsF(f,:);
        end
        h  = scatter(ax, xx, yy, mSize, ...
            'Marker',          markersF{f}, ...
            'MarkerFaceColor', dotColor, ...
            'MarkerEdgeColor', dotColor, ...
            'MarkerFaceAlpha', alphaVal, ...
            'MarkerEdgeAlpha', alphaVal, ...
            'DisplayName',     char(folderCats(f)));
        if ~hasHandle(f), lgHandles(f) = h; hasHandle(f) = true; end
    end
    m = mean(yAll, 'omitnan');
    plot(ax, [xPos(k)-meanHalf, xPos(k)+meanHalf], [m, m], 'k-', 'LineWidth', 1.8);
end

% =========================================================================
% 12. Fix axes and stretch background bands
% =========================================================================
margin = xSpacing * 0.5;
xlim(ax, [xPos(1)-margin, xPos(end)+margin]);
set(fig, 'Renderer', 'opengl');

validMetric = metricF(~isnan(metricF));
padAll = 0.08 * range(validMetric);
if ~isfinite(padAll) || padAll == 0
    padAll = 0.05 * max(abs(validMetric));
    if ~isfinite(padAll), padAll = 1; end
end
yl = ylim(ax);
if ~isempty(YLimMaxSet), yl(2) = YLimMaxSet; end
if ~isempty(YLimMinSet), yl(1) = YLimMinSet; end
ylim(ax, [yl(1), yl(2) + padAll]);
yl = ylim(ax);

% % Stretch bands to full y-extent
% for k = 1:K
%     b = findobj(ax, 'Tag', sprintf('band_%d', k));
%     if ~isempty(b)
%         set(b, 'YData', [yl(1) yl(1) yl(2) yl(2)]);
%     end
% end

% Send bands and violins to the back
% bands = findobj(ax, '-regexp', 'Tag', '^band_');
% if ~isempty(bands), uistack(bands, 'bottom'); end
if showViolin
    violinPatches = findobj(ax, '-regexp', 'Tag', '^violin');
    if ~isempty(violinPatches), uistack(violinPatches, 'bottom'); end
end

% =========================================================================
% 13. Compact significance brackets (black, significant pairs only)
%     Algorithm:
%       1. Filter to significant pairs whose both endpoints are in K groups.
%       2. Sort by span width (ascending) so short brackets are drawn first
%          and sit lower.
%       3. For each bracket, place it just above the tallest bracket top
%          already occupying any column in its span.
% =========================================================================
if ~isempty(sigPairs)
    % Map "<compare>_<level>" -> x position
    % compare1 groups sit at k-dx, compare2 groups at k+dx
    groupXPos = containers.Map('KeyType','char','ValueType','double');
    groupIdxMap = containers.Map('KeyType','char','ValueType','double');
    for k = 1:K
        groupXPos(char(levelsPresent{k}))   = xPos(k);  % real x-position
        groupIdxMap(char(levelsPresent{k})) = k;         % integer index for colTop
    end

pairInfo = [];
    for r = 1:numel(sigPairs)
        gA = sigPairs{r}{1};  gB = sigPairs{r}{2};
        if ~isKey(groupXPos, gA) || ~isKey(groupXPos, gB)
            fprintf('WARNING: skipped pair %s vs %s — not in groupXPos\n', gA, gB);
            continue;
        end
        xA = groupXPos(gA);   xB = groupXPos(gB);
        pairInfo(end+1).xA   = min(xA, xB);
        pairInfo(end).xB     = max(xA, xB);
        pairInfo(end).iA     = min(groupIdxMap(gA), groupIdxMap(gB));
        pairInfo(end).iB     = max(groupIdxMap(gA), groupIdxMap(gB));
        pairInfo(end).star   = sigPairs{r}{3};
        pairInfo(end).span   = abs(xB - xA);
    end

    if ~isempty(pairInfo)
        [~, sortOrd] = sort([pairInfo.span]);
        pairInfo = pairInfo(sortOrd);

        % Initialise per-column ceiling using actual data max per group
        colTop = zeros(1, K);
        for k = 1:K
            grpVals = metricF(groupF == levelsPresent{k});
            colTop(k) = max(grpVals(~isnan(grpVals)), [], 'omitnan');
            if ~isfinite(colTop(k))
                yl = ylim(ax);
                colTop(k) = yl(1);
            end
        end
        % add violin clearance
        yl = ylim(ax);
        colTop = colTop + 0.04 * (yl(2) - yl(1));

        bracketStep  = 0.07 * (yl(2) - yl(1));
        bracketColor = [0 0 0];

        for r = 1:numel(pairInfo)
            xA   = pairInfo(r).xA;
            xB   = pairInfo(r).xB;
            star = pairInfo(r).star;
            iA     = pairInfo(r).iA;
            iB     = pairInfo(r).iB;

            % columns spanned: xA and xB are integer indices now
            inSpan = iA:iB;
            yBase  = max(colTop(inSpan)) + bracketStep;

            yl = ylim(ax);
            if yBase > yl(2)
                ylim(ax, [yl(1), yBase + bracketStep]);
                yl = ylim(ax);
            end

            footGap = 0.025 * (yl(2) - yl(1));
            yFootA  = colTop(iA) + footGap;
            yFootB  = colTop(iB) + footGap;
            plot(ax, [xA xA], [yFootA, yBase], 'Color', bracketColor, 'LineWidth', 1.2);
            plot(ax, [xB xB], [yFootB, yBase], 'Color', bracketColor, 'LineWidth', 1.2);
            plot(ax, [xA xB], [yBase,  yBase], 'Color', bracketColor, 'LineWidth', 1.2);
            text(ax, (xA+xB)/2, yBase + 0.005*(yl(2)-yl(1)), star, ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', bracketColor);

            % update colTop for all columns spanned by this bracket
            colTop(inSpan) = yBase;
        end

        yl = ylim(ax);
        ylim(ax, [yl(1), yl(2) + 0.04*(yl(2)-yl(1))]);
    end
end

% if ~isempty(sigPairs)
%     % Map group name -> x index
%     groupXMap = containers.Map(groupsPresent, num2cell(1:K));
% 
%     % Build valid pair list with x indices and span width
%     pairInfo = [];
%     for r = 1:numel(sigPairs)
%         gA = sigPairs{r}{1};  gB = sigPairs{r}{2};
%         if ~isKey(groupXMap, gA) || ~isKey(groupXMap, gB), continue; end
%         xA = groupXMap(gA);   xB = groupXMap(gB);
%         pairInfo(end+1).xA   = min(xA, xB);           %#ok<AGROW>
%         pairInfo(end).xB     = max(xA, xB);
%         pairInfo(end).star   = sigPairs{r}{3};
%         pairInfo(end).span   = abs(xB - xA);
%     end
% 
%     if ~isempty(pairInfo)
%         % Sort by span ascending (short brackets drawn lower)
%         [~, sortOrd] = sort([pairInfo.span]);
%         pairInfo = pairInfo(sortOrd);
% 
%         % Initialise per-column ceiling (highest y used above each tick)
%         colTop = zeros(1, K);
%         for k = 1:K
%             grpVals = metricF(groupF == groupsPresent{k});
%             colTop(k) = max(grpVals(~isnan(grpVals)), [], 'omitnan');
%             if ~isfinite(colTop(k))
%                 yl = ylim(ax); 
%                 colTop(k) = yl(2);
%             end
%         end
% 
%         yl = ylim(ax);
%         bracketStep = 0.07 * (yl(2) - yl(1));
%         bracketColor = [0 0 0];
% 
%         for r = 1:numel(pairInfo)
%             xA   = pairInfo(r).xA;
%             xB   = pairInfo(r).xB;
%             star = pairInfo(r).star;
% 
%             % y base = one step above the tallest column in the span
%             yBase = max(colTop(xA:xB)) + bracketStep;
% 
%             % Expand axes if needed
%             yl = ylim(ax);
%             if yBase > yl(2)
%                 ylim(ax, [yl(1), yBase + bracketStep]);
%                 yl = ylim(ax);
%             end
% 
%             tickH = 0.012 * (yl(2) - yl(1));
% 
%             % Draw bracket: two vertical ticks + horizontal bar
%             plot(ax, [xA xA], [yBase - tickH, yBase], ...
%                 'Color', bracketColor, 'LineWidth', 1.2);
%             plot(ax, [xB xB], [yBase - tickH, yBase], ...
%                 'Color', bracketColor, 'LineWidth', 1.2);
%             plot(ax, [xA xB],  [yBase, yBase],          ...
%                 'Color', bracketColor, 'LineWidth', 1.2);
%             text(ax, (xA + xB) / 2, yBase + 0.005 * (yl(2) - yl(1)), star, ...
%                 'HorizontalAlignment', 'center', ...
%                 'FontSize', 10, 'Color', bracketColor);
% 
%             % Update column ceilings over the entire span
%             colTop(xA:xB) = yBase;
%         end
% 
%         % Small extra padding above the highest bracket
%         yl = ylim(ax);
%         ylim(ax, [yl(1), yl(2) + 0.04 * (yl(2) - yl(1))]);
%     end
% end

% =========================================================================
% 14. N labels
% =========================================================================
yl = ylim(ax);
nLabelPad = 0.01 * (yl(2) - yl(1));

for k = 1:K
    yl   = ylim(ax);
    yTop = yl(2) + nLabelPad;
    if nPerGroup(k) > 0
        text(ax, xPos(k), yTop, sprintf('n=%d', nPerGroup(k)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'bottom', ...
            'FontSize', 7, 'Color', [0.3 0.3 0.3], 'Clipping', 'off');
    end
end

yl = ylim(ax);
ylim(ax, [yl(1), yl(2) + 0.05 * (yl(2) - yl(1))]);

% =========================================================================
% 15. Axes decoration
% =========================================================================
xticks(ax, xPos);

% Shorten x-axis labels to at most 3 underscore-separated tokens
shortLabels = cellfun(@shortenLabel3, levelsPresent, 'UniformOutput', false);
xticklabels(ax, shortLabels);
ax.XTickLabelRotation = 30;
ax.TickLabelInterpreter = 'none';

ylabel(ax, sprintf('%s %s', varNameStr, unitName), 'Interpreter', 'none');
title(ax, varNameStr, 'Interpreter', 'none');

subtitleParts = {sprintf('one-way ANOVA, post-hoc: %s', postHocTest)};
if isfinite(CutoffBottom), subtitleParts{end+1} = sprintf('cut_bot=%.4g', CutoffBottom); end
if isfinite(CutoffTop),    subtitleParts{end+1} = sprintf('cut_top=%.4g', CutoffTop);    end
subtitle(ax, strjoin(subtitleParts, '  |  '), 'FontSize', 8, 'Color', [0.4 0.4 0.4]);

grid(ax, 'on');
set(ax, 'Box', 'off', 'Layer', 'top');

% Legend (DataFolder, shortened)
useIdx = find(hasHandle);
if ~isempty(useIdx)
    labels      = cellstr(folderCats(useIdx));
    labelsShort = cellfun(@shortenLabel3, labels, 'UniformOutput', false);
    legend(lgHandles(useIdx), labelsShort, ...
        'Location', 'southoutside', 'Interpreter', 'none', ...
        'NumColumns', min(nF, 3));
end

% =========================================================================
% 16. Overview table
% =========================================================================
overviewTbl = table( ...
    string(levelsPresent(:)), ...
    zeros(K,1), zeros(K,1), zeros(K,1), ...
    'VariableNames', {'Group','Mean','SD','N'});

for k = 1:K
    grp  = levelsPresent{k};
    vals = metricF(groupF == grp);
    vals = vals(~isnan(vals));
    overviewTbl.Mean(k) = mean(vals, 'omitnan');
    overviewTbl.SD(k)   = std(vals,  'omitnan');
    overviewTbl.N(k)    = numel(vals);
end

fprintf('\n--- Overview of %s per group ---\n', varNameStr);
disp(overviewTbl);

hold(ax, 'off');
end


% =========================================================================
% LOCAL HELPERS
% =========================================================================

function out = pToStars(p)
    if p < 0.001,     out = '***';
    elseif p < 0.01,  out = '**';
    elseif p < 0.05,  out = '*';
    else,             out = 'ns';
    end
end

function out = shortenLabel3(s)
    parts = strings(1,3);
    rest  = char(s);
    for i = 1:3
        [tok, rest] = strtok(rest, '_');
        if isempty(tok)
            parts = parts(1:i-1);
            break
        end
        parts(i) = string(tok);
        if isempty(rest)
            parts = parts(1:i);
            break
        end
        rest = rest(2:end);
    end
    out = char(strjoin(parts, '_'));
end

function drawViolin(ax, data, xCenter, halfWidth, nKde, faceColor, faceAlpha, tagStr)
    if numel(data) < 3, return; end
    dataMin = min(data);  dataMax = max(data);
    span    = dataMax - dataMin;
    if span == 0
        yq  = dataMin + [-1e-9; 1e-9];
        den = [1; 1];
    else
        pad = 0.10 * span;
        yq  = linspace(dataMin - pad, dataMax + pad, nKde)';
        den = ksdensity(data, yq);
    end
    maxDen = max(den);
    if maxDen == 0, return; end
    denNorm = (den / maxDen) * halfWidth;
    xRight  = xCenter + denNorm;
    xLeft   = xCenter - denNorm;
    xPoly   = [xRight; flipud(xLeft)];
    yPoly   = [yq;     flipud(yq)];
    patch(ax, xPoly, yPoly, faceColor, ...
        'FaceAlpha', faceAlpha, ...
        'EdgeColor', faceColor * 0.6, ...
        'EdgeAlpha', min(faceAlpha + 0.3, 1), ...
        'LineWidth', 0.8, ...
        'Tag',       tagStr);
end