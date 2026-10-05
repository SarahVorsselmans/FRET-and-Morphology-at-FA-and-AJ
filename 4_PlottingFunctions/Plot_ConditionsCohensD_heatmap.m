function fig = Plot_ConditionsCohensD_heatmap(comparisonTable, metrics, options)
% PLOT_COHENSD_HEATMAP
%   Visualizes Cohen's D from Calc_compareConditionsToReference output
%   as a diverging heatmap (conditions × metrics).
%
% INPUTS:
%   comparisonTable - output of Calc_compareConditionsToReference
%   metrics         - same cell array of metric names used to build the table
%
% OPTIONS:
%   CLim            - [min max] color axis limits. Default: symmetric around max |d|
%   Annotate        - show numeric values in cells (default: true)
%   FontSize        - cell annotation font size (default: 9)
%   Title           - figure title string (default: "Cohen's D heatmap")
%   ColorMap        - Nx3 colormap matrix (default: built-in diverging RdBu)

arguments
    comparisonTable     table
    metrics             cell
    options.Levels      = []
    options.CLim        (1,2) double = [NaN NaN]
    options.Annotate    logical      = true
    options.FontSize    double       = 9
    options.Title       string       = "Cohen's D — Conditions vs. Reference"
    options.ColorMap    double       = []
    options.TitleInfo   string       = ""
end

% --- 0. Select / reorder condition rows ---
if isempty(options.Levels)
    rowOrder = 1:height(comparisonTable);
else
    levels = cellstr(options.Levels);
    allLabels = cellstr(string(comparisonTable.conditionLabel));
    rowOrder = NaN(numel(levels), 1);
    for i = 1:numel(levels)
        idx = find(strcmp(allLabels, levels{i}));
        if isempty(idx)
            error('Plot_ConditionsCohensD_heatmap:LevelNotFound', ...
                'conditionLabel "%s" not found in comparisonTable.', levels{i});
        elseif numel(idx) > 1
            error('Plot_ConditionsCohensD_heatmap:LevelNotUnique', ...
                'conditionLabel "%s" matches multiple rows in comparisonTable.', levels{i});
        end
        rowOrder(i) = idx;
    end
end
comparisonTable = comparisonTable(rowOrder, :);

% --- 1. Extract Cohen's D matrix [nConds x nMetrics] ---
nMetrics = numel(metrics);
nConds   = height(comparisonTable);

D = NaN(nConds, nMetrics);
for m = 1:nMetrics
    colName  = sprintf('%s_cohensD', metrics{m});
    D(:, m)  = comparisonTable.(colName);
end

condLabels   = comparisonTable.conditionLabel;   % cell or string array
metricLabels = metrics;

% --- 2. Color limits ---
maxAbsD = max(abs(D(:)), [], 'omitnan');
if any(isnan(options.CLim))
    climRange = [-maxAbsD, maxAbsD];
else
    climRange = options.CLim;
end

% --- 3. Colormap: diverging blue–white–red ---
if isempty(options.ColorMap)
    cmap = divergingColormap(256);
else
    cmap = options.ColorMap;
end

% --- 4. Draw heatmap ---
fig = figure('Color', 'w');
ax  = axes(fig);
imagesc(ax, D);
colormap(ax, cmap);
clim(ax, climRange);
cb       = colorbar(ax);
cb.Label.String = "Cohen's d";

% Axes labels
ax.XTick      = 1:nMetrics;
ax.XTickLabel = metricLabels;
ax.XTickLabelRotation = 30;
ax.YTick      = 1:nConds;
ax.YTickLabel = condLabels;
ax.TickLength = [0 0];
ax.FontSize   = 10;
ax.TickLabelInterpreter = 'none';
axis(ax, 'tight');

title(ax, sprintf('%s %s', options.Title, options.TitleInfo), 'FontSize', 12, 'FontWeight', 'bold');

% --- 5. Cell annotations ---
if options.Annotate
    for r = 1:nConds
        for c = 1:nMetrics
            val = D(r, c);
            if ~isnan(val)
                % Dark text on light cells, light text on dark cells
                normVal  = (val - climRange(1)) / (climRange(2) - climRange(1));  % 0–1
                txtColor = pickTextColor(normVal);
                text(ax, c, r, sprintf('%.2f', val), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment',   'middle', ...
                    'FontSize',            options.FontSize, ...
                    'Color',               txtColor);
            end
        end
    end
end

% --- 6. Grid lines between cells ---
hold(ax, 'on');
for c = 1.5:nMetrics - 0.5
    xline(ax, c, 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
end
for r = 1.5:nConds - 0.5
    yline(ax, r, 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
end
hold(ax, 'off');
end


%% ---- Local helpers ------------------------------------------------
function cmap = divergingColormap(n)
    % Blue–white–red diverging map
    half  = floor(n/2);
    blue  = [linspace(0.02, 1, half)', linspace(0.18, 1, half)', linspace(0.47, 1, half)'];
    red   = [linspace(1, 0.70, half)', linspace(1, 0.09, half)', linspace(1, 0.09, half)'];
    cmap  = [blue; red];
end

function c = pickTextColor(normVal)
    % White text when cell is very dark (near 0 or 1), black otherwise
    if normVal < 0.25 || normVal > 0.75
        c = [1 1 1];   % white
    else
        c = [0 0 0];   % black
    end
end

%% More difference between colors (more white)
% function fig = Plot_ConditionsCohensD_heatmap(comparisonTable, metrics, options)
% % PLOT_COHENSD_HEATMAP
% %   Visualizes Cohen's D from Calc_compareConditionsToReference output
% %   as a diverging heatmap (conditions × metrics).
% %
% % INPUTS:
% %   comparisonTable - output of Calc_compareConditionsToReference
% %   metrics         - same cell array of metric names used to build the table
% %
% % OPTIONS:
% %   CLim            - [min max] color axis limits. Default: symmetric around max |d|
% %   Annotate        - show numeric values in cells (default: true)
% %   FontSize        - cell annotation font size (default: 9)
% %   Title           - figure title string (default: "Cohen's D heatmap")
% %   ColorMap        - Nx3 colormap matrix (default: built-in diverging RdBu)
% 
% arguments
%     comparisonTable     table
%     metrics             cell
%     options.CLim        (1,2) double = [NaN NaN]
%     options.Annotate    logical      = true
%     options.FontSize    double       = 9
%     options.Title       string       = "Cohen's D — Conditions vs. Reference"
%     options.ColorMap    double       = []
%     options.TitleInfo   string       = ""
% end
% 
% % --- 1. Extract Cohen's D matrix [nConds x nMetrics] ---
% nMetrics = numel(metrics);
% nConds   = height(comparisonTable);
% 
% D = NaN(nConds, nMetrics);
% for m = 1:nMetrics
%     colName  = sprintf('%s_cohensD', metrics{m});
%     D(:, m)  = comparisonTable.(colName);
% end
% 
% condLabels   = comparisonTable.conditionLabel;   % cell or string array
% metricLabels = metrics;
% 
% % --- 2. Color limits ---
% maxAbsD = max(abs(D(:)), [], 'omitnan');
% if any(isnan(options.CLim))
%     climRange = [-maxAbsD, maxAbsD];
% else
%     climRange = options.CLim;
% end
% 
% % --- 3. Colormap: diverging blue–white–red ---
% if isempty(options.ColorMap)
%     cmap = divergingColormap(256, climRange, 0.5); % colors only very visible from absolute value of 0.5
% else
%     cmap = options.ColorMap;
% end
% 
% % --- 4. Draw heatmap ---
% fig = figure('Color', 'w');
% ax  = axes(fig);
% 
% imagesc(ax, D);
% colormap(ax, cmap);
% clim(ax, climRange);
% cb       = colorbar(ax);
% cb.Label.String = "Cohen's d";
% 
% % Axes labels
% ax.XTick      = 1:nMetrics;
% ax.XTickLabel = metricLabels;
% ax.XTickLabelRotation = 30;
% ax.YTick      = 1:nConds;
% ax.YTickLabel = condLabels;
% ax.TickLength = [0 0];
% ax.FontSize   = 10;
% ax.TickLabelInterpreter = 'none';
% axis(ax, 'tight');
% 
% title(ax, sprintf('%s %s', options.Title, options.TitleInfo), 'FontSize', 12, 'FontWeight', 'bold');
% 
% % --- 5. Cell annotations ---
% if options.Annotate
%     for r = 1:nConds
%         for c = 1:nMetrics
%             val = D(r, c);
%             if ~isnan(val)
%                 % Dark text on light cells, light text on dark cells
%                 normVal  = (val - climRange(1)) / (climRange(2) - climRange(1));  % 0–1
%                 txtColor = pickTextColor(normVal);
%                 text(ax, c, r, sprintf('%.2f', val), ...
%                     'HorizontalAlignment', 'center', ...
%                     'VerticalAlignment',   'middle', ...
%                     'FontSize',            options.FontSize, ...
%                     'Color',               txtColor);
%             end
%         end
%     end
% end
% 
% % --- 6. Grid lines between cells ---
% hold(ax, 'on');
% for c = 1.5:nMetrics - 0.5
%     xline(ax, c, 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
% end
% for r = 1.5:nConds - 0.5
%     yline(ax, r, 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
% end
% hold(ax, 'off');
% end
% 
% 
% %% ---- Local helpers ------------------------------------------------
% 
% function cmap = divergingColormap(n, climRange, threshold)
% % Blue–white–red diverging map with a white zone for |value| < threshold)
% 
% if nargin < 3
%     threshold = 0; % default: no dead zone
% end
% 
% % Fraction of the full color axis range that should stay white
% range = climRange(2) - climRange(1);
% whiteFrac = min(1, (2 * threshold) / range);  % proportion of [climRange] within +/-threshold
% 
% % Colormap
% half  = floor(n/2);
% % blue  = [linspace(0.02, 1, half)', linspace(0.18, 1, half)', linspace(0.47, 1, half)'];
% % red   = [linspace(1, 0.70, half)', linspace(1, 0.09, half)', linspace(1, 0.09, half)'];
% % cmap  = [blue; red];
% nWhiteHalf = round((whiteFrac/2) * n);   % white samples per side
% nColorHalf = half - nWhiteHalf;
% 
% if nColorHalf < 1
%     nColorHalf = 1;
% end
% 
% % Blue side: blue -> white
% blue = [linspace(0.02, 1, nColorHalf)', linspace(0.18, 1, nColorHalf)', linspace(0.47, 1, nColorHalf)'];
% % White padding (so the inner band stays exactly white)
% blueWhitePad = repmat([1 1 1], nWhiteHalf, 1);
% 
% % Red side: white -> red
% red = [linspace(1, 0.70, nColorHalf)', linspace(1, 0.09, nColorHalf)', linspace(1, 0.09, nColorHalf)'];
% redWhitePad = repmat([1 1 1], nWhiteHalf, 1);
% 
% cmap = [blue; blueWhitePad; redWhitePad; red];
% 
% % Pad/truncate to exactly n rows in case of rounding
% if size(cmap,1) < n
%     cmap = [cmap; repmat(cmap(end,:), n - size(cmap,1), 1)];
% elseif size(cmap,1) > n
%     cmap = cmap(1:n, :);
% end
% end
% 
% % function cmap = divergingColormap(n)
% % % Blue–white–red diverging map
% %     half  = floor(n/2);
% %     blue  = [linspace(0.02, 1, half)', linspace(0.18, 1, half)', linspace(0.47, 1, half)'];
% %     red   = [linspace(1, 0.70, half)', linspace(1, 0.09, half)', linspace(1, 0.09, half)'];
% %     cmap  = [blue; red];
% % end
% 
% %%
% function c = pickTextColor(normVal)
% % White text when cell is very dark (near 0 or 1), black otherwise
% if normVal < 0.25 || normVal > 0.75
%     c = [1 1 1];   % white
% else
%     c = [0 0 0];   % black
% end
% end