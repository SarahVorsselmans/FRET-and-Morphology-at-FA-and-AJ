function fig = Plot_CompareVar(tableT, XVar, YVar, options)
% could be: [fig, axs, stats] = Plot_CompareVar(...)
%PLOT_COMPAREVAR  Compare two variables across condition combinations.
%
% [fig, axs, stats] = Plot_CompareVar(tableT, XVar, YVar, options)
%
% ── INPUTS ────────────────────────────────────────────────────────────────
%   tableT   : table where each row = one ROI.
%              Each variable column holds a Nx1 double (per-object values).
%              Must contain categorical columns for Cond1 and Cond2.
%
%   XVar     : string — column name for X-axis variable
%   YVar     : string — column name for Y-axis variable
%
%   options  : name-value pairs (all optional, sensible defaults provided)
%
% ── OPTIONS ───────────────────────────────────────────────────────────────
%   Cond1        : string, column name for disease category  (def: 'diseaseCat')
%   Cond2        : column name for extra-info category       (def: 'extrainfoCat')
%   ExpCat       : column name for biological replicates     (def: 'expCat')
%   Levels       : cell array of Cond2 levels to include     (def: all levels)
%
%   PlotType     : 'scatter' | 'binned'                      (def: 'scatter')
%   NBins        : number of X bins for binned plot          (def: 10)
%
%   ScatterMode  : 'roi' | 'all'                             (def: 'roi')
%                  'roi'  → one dot per ROI (per-ROI mean of Nx1 vector)
%                  'all'  → one dot per object (all raw values in Nx1)
%
%   ShowDots     : logical — show individual dots in scatter plot (def: true)
%   ShowMeans    : logical — overlay per-expCat mean markers     (def: true)
%   ShowEllipse  : logical — show 5%-95% percentile ellipse per
%                  condition in scatter plot                 (def: true)
%   ScatterStyle : string = "points" | "heatmap"             (def: heatmap)
%   EllipseStats : logical — compare ellipses across conditions  (def: false)
%                  (only relevant in combined / cross-condition panels)
%
%   Stats        : logical — compute & display statistics    (def: true)
%   StatTest     : 'spearman' | 'pearson'                    (def: 'spearman')
%                  -> rho: How strongly X and Y change together 
%                  (+1 = pos relation, 0 = no relation, -1 = neg relation)
%                  -> p-value: Whether the observed correlation is likely real or just noise
%                  (< 0.05 = significant, > 0.05 = not siginficant)
%
%   CondColors   : struct with fields = Cond1 values, each a struct with
%                  fields = Cond2 values, each a palette name string.
%                  Example:
%                    options.CondColors.CT.rest  = 'blues';
%                    options.CondColors.CT.R1    = 'indigos';
%                    options.CondColors.CCM.rest = 'redes';
%                  If empty, auto-assigns distinct palettes.
%
%   FigSize      : [width height] in pixels                  (def: auto)
%   AxFontSize   : font size for axes                        (def: 9)
%   MarkerSize   : marker size for scatter dots              (def: 18)
%   MeanMarkerSize: marker size for mean markers             (def: 60)
%   Alpha        : transparency for scatter dots             (def: 0.4)
%
% ── OUTPUTS ───────────────────────────────────────────────────────────────
%   fig   : figure handle
%   axs   : struct with fields matching panel labels, each an axes handle
%   stats : struct with fields matching panel labels, each containing:
%             .rho  (Spearman or Pearson correlation)
%             .pval
%             .n    (number of ROIs / pooled objects)
%             .ellipseTest (if EllipseStats=true and ≥2 conditions)
%
% ── PANEL LAYOUT ──────────────────────────────────────────────────────────
%   Given Cond1 = {A, B} and Levels = {L1, L2, L3}:
%
%   Row 1: A_L1 | A_L2 | A_L3 | A_L1+A_L2+A_L3  (per expCat for binned)
%   Row 2: B_L1 | B_L2 | B_L3 | B_L1+B_L2+B_L3  (per expCat for binned)
%   Row 3: A_L1+B_L1 | A_L2+B_L2 | A_L3+B_L3 | all  (combined for binned)
%
% ── DATA NOTES ────────────────────────────────────────────────────────────
%   Scatter 'roi'  mode: one dot per ROI = mean of that ROI's Nx1 vector.
%   Scatter 'all'  mode: one dot per object = every value in each Nx1 vector.
%   Binned  ribbon: always uses all raw values (pooled Nx1 across ROIs),
%                   with per-ROI equal weighting so no single ROI dominates.

arguments
    tableT                      table
    XVar                        string
    YVar                        string
    options.Cond1               string   = "diseaseCat"
    options.Cond2               string   = "extrainfoCat"
    options.ExpCat              string   = "expCat"
    options.Levels              cell     = {}
    options.PlotType            string   = "scatter"
    options.NBins               double   = 10
    options.ScatterMode         string   = "roi"
    options.ShowDots            logical  = true
    options.ShowMeans           logical  = true
    options.ShowEllipse         logical  = true
    options.ScatterStyle        string   = "heatmap"
    options.globalYLim          double   = []
    options.EllipseStats        logical  = false
    options.Stats               logical  = true
    options.StatTest            string   = "spearman"
    options.CondColors                   = []
    options.FigSize             double   = []
    options.AxFontSize          double   = 9
    options.MarkerSize          double   = 18
    options.MeanMarkerSize      double   = 60
    options.Alpha               double   = 0.4
end

% ══════════════════════════════════════════════════════════════════════════
%  0. VALIDATE & RESOLVE LEVELS
% ══════════════════════════════════════════════════════════════════════════
cond1Col = options.Cond1;
cond2Col = options.Cond2;
expCol   = options.ExpCat;

% Validate PlotType
if ~ismember(lower(options.PlotType), {'scatter','binned'})
    error('Plot_CompareVar: PlotType must be ''scatter'' or ''binned''.');
end
% Validate ScatterMode
if ~ismember(lower(options.ScatterMode), {'roi','all'})
    error('Plot_CompareVar: ScatterMode must be ''roi'' or ''all''.');
end

% Get unique Cond1 values (e.g. CT, CCM)
c1Vals = unique(string(tableT.(cond1Col)), 'stable');
nC1    = numel(c1Vals);

% Get unique Cond2 levels (e.g. rest, R1, R2)
if isempty(options.Levels)
    allLevels = unique(string(tableT.(cond2Col)), 'stable');
else
    allLevels = string(options.Levels(:)');
end
nLev = numel(allLevels);

% Number of columns = nLev individual + 1 combined
nCols = nLev + 1;
% Number of rows    = nC1 individual + 1 cross-condition
nRows = nC1 + 1;

% ══════════════════════════════════════════════════════════════════════════
%  1. BUILD PANEL DEFINITIONS
%     Each panel = struct with .label, .conditions, .isCombined
% ══════════════════════════════════════════════════════════════════════════
panels = buildPanelDefs(c1Vals, allLevels, nRows, nCols);

% ══════════════════════════════════════════════════════════════════════════
%  2. RESOLVE COLORS per (Cond1, Cond2) combination
% ══════════════════════════════════════════════════════════════════════════
colorMap = resolveColorMap(c1Vals, allLevels, options.CondColors);

% ══════════════════════════════════════════════════════════════════════════
%  3. CREATE FIGURE
% ══════════════════════════════════════════════════════════════════════════
panelW = 200; panelH = 190;
padL = 55; padR = 20; padT = 40; padB = 50;
gapX = 45; gapY = 55;

if isempty(options.FigSize)
    figW = padL + nCols*panelW + (nCols-1)*gapX + padR;
    figH = padT + nRows*panelH + (nRows-1)*gapY + padB;
else
    figW = options.FigSize(1);
    figH = options.FigSize(2);
end

fig = figure('Color','w', 'Units','pixels', 'Position',[100 100 figW figH]);

% ══════════════════════════════════════════════════════════════════════════
%  4. FIRST PASS — collect all condData and compute global axis limits
% ══════════════════════════════════════════════════════════════════════════
allCondData = cell(nRows, nCols);
for r = 1:nRows
    for c = 1:nCols
        pDef = panels{r,c};
        if isempty(pDef), continue; end
        allCondData{r,c} = collectConditionData(tableT, XVar, YVar, ...
                               pDef.conditions, cond1Col, cond2Col, expCol);
    end
end

% Derive limits
globalXLim = [Inf, -Inf];
globalYLim = [Inf, -Inf];

switch lower(options.PlotType)
    case 'scatter'
        % Use raw values from the 'all' panel (bottom-right)
        allPanelData = allCondData{nRows, nCols};
        if ~isempty(allPanelData)
            for k = 1:numel(allPanelData)
                cd = allPanelData(k);
                if strcmpi(options.ScatterMode, 'all')
                    xVals = cd.poolX;  yVals = cd.poolY;
                else
                    xVals = cd.roiX;   yVals = cd.roiY;
                end
                xVals = xVals(isfinite(xVals));
                yVals = yVals(isfinite(yVals));
                if ~isempty(xVals)
                    globalXLim(1) = min(globalXLim(1), min(xVals));
                    globalXLim(2) = max(globalXLim(2), max(xVals));
                end
                if ~isempty(yVals)
                    globalYLim(1) = min(globalYLim(1), min(yVals));
                    globalYLim(2) = max(globalYLim(2), max(yVals));
                end
            end
        end

    case 'binned'
        % Step 1: get xLim from raw data to define bin edges
        allPanelData = allCondData{nRows, nCols};
        if ~isempty(allPanelData)
            for k = 1:numel(allPanelData)
                xVals = allPanelData(k).poolX;
                xVals = xVals(isfinite(xVals));
                if ~isempty(xVals)
                    globalXLim(1) = min(globalXLim(1), min(xVals));
                    globalXLim(2) = max(globalXLim(2), max(xVals));
                end
            end
        end
        xPad = max((globalXLim(2)-globalXLim(1))*0.05, eps);
        globalXLim = globalXLim + [-xPad, xPad];

        % Step 2: compute bin edges, then derive both xLim and yLim
        %         from the actual binned statistics (mean ± SEM).
        %         This avoids raw outliers inflating the axes.
        edges    = linspace(globalXLim(1), globalXLim(2), options.NBins+1);
        binCtrs  = (edges(1:end-1) + edges(2:end)) / 2;
        xOccupied = [Inf, -Inf];   % tightest x range that has valid bins

        for r = 1:nRows
            for c = 1:nCols
                condData = allCondData{r,c};
                if isempty(condData), continue; end
                for k = 1:numel(condData)
                    cd = condData(k);
                    valid = isfinite(cd.poolX) & isfinite(cd.poolY) & isfinite(cd.poolW);
                    if sum(valid) < 2, continue; end
                    [mn, se] = weightedBinStats( ...
                        cd.poolX(valid), cd.poolY(valid), cd.poolW(valid), edges);

                    validBins = isfinite(mn) & isfinite(se);
                    if ~any(validBins), continue; end

                    % y: span of mean ± SEM
                    yLow  = mn(validBins) - se(validBins);
                    yHigh = mn(validBins) + se(validBins);
                    globalYLim(1) = min(globalYLim(1), min(yLow));
                    globalYLim(2) = max(globalYLim(2), max(yHigh));

                    % x: only bins that actually contain data
                    xOccupied(1) = min(xOccupied(1), min(binCtrs(validBins)));
                    xOccupied(2) = max(xOccupied(2), max(binCtrs(validBins)));
                end
            end
        end

        % Replace globalXLim with the tighter occupied range
        if isfinite(xOccupied(1))
            globalXLim = xOccupied;
        end
end

% Add 5% padding
xPad = max((globalXLim(2)-globalXLim(1))*0.05, eps);
yPad = max((globalYLim(2)-globalYLim(1))*0.05, eps);
if isfinite(globalXLim(1)), globalXLim = globalXLim + [-xPad, xPad]; end
if isfinite(globalYLim(1)), globalYLim = globalYLim + [-yPad, yPad]; end

% ══════════════════════════════════════════════════════════════════════════
%  5. SECOND PASS — create axes, plot, annotate, apply uniform limits
% ══════════════════════════════════════════════════════════════════════════
axs   = struct();
stats = struct();

for r = 1:nRows
    for c = 1:nCols
        pDef     = panels{r,c};
        condData = allCondData{r,c};
        if isempty(pDef) || isempty(condData), continue; end

        % isCombined: last column of top rows OR entire last row
        isCombined = pDef.isCombined;

        % Create axes
        xPos = (padL + (c-1)*(panelW+gapX)) / figW;
        yPos = (padB + (nRows-r)*(panelH+gapY)) / figH;
        wPos = panelW / figW;
        hPos = panelH / figH;
        ax = axes('Parent', fig, 'Units','normalized', ...
                  'Position',[xPos yPos wPos hPos], ...
                  'FontSize', options.AxFontSize, ...
                  'Box','on', 'TickDir','out');
        hold(ax,'on');

        % Plot
        switch lower(options.PlotType)
            case 'scatter'
                ellSt = plotScatter(ax, condData, colorMap, options, isCombined);
            case 'binned'
                ellSt = [];
                ax    = plotBinned(ax, condData, colorMap, options, ...
                                   globalXLim, isCombined);
        end

        % Apply uniform axis limits
        if isfinite(globalXLim(1)), xlim(ax, globalXLim); end
        if isfinite(globalYLim(1)), ylim(ax, globalYLim); end

        % Correlation statistics (not for combined panels to avoid repetition)
        st = [];
        isCombinedCol = (c == nCols);
        isCombinedRow = (r == nRows);
        showStatsHere = options.Stats && ~(isCombinedCol || isCombinedRow);
        if showStatsHere
            st = computeStats(condData, options.StatTest, ...
                              options.PlotType, options.ScatterMode,...
                              globalXLim, options.NBins);
            addStatsAnnotation(ax, st, options.AxFontSize);
            % also compute binned slope stats for binned plots
            if strcmpi(options.PlotType, 'binned')
                slopeSt = computeBinnedSlope(condData, globalXLim, options.NBins);
                addSlopeStatsAnnotation(ax, slopeSt, options.AxFontSize);
                % merge into st so it's returned in the output struct too
                fn = fieldnames(slopeSt);
                for fi = 1:numel(fn)
                    st.(fn{fi}).slope    = slopeSt.(fn{fi}).slope;
                    st.(fn{fi}).slope_se = slopeSt.(fn{fi}).slope_se;
                    st.(fn{fi}).slope_p  = slopeSt.(fn{fi}).pval;
                end
            end
        end

        % Attach ellipse comparison stats if computed
        if ~isempty(ellSt)
            if isempty(st), st = struct(); end
            st.ellipseTest = ellSt;
        end

        % Labels
        if r == nRows
            xlabel(ax, strrep(XVar,'_',' '), 'FontSize', options.AxFontSize);
        end
        if c == 1
            ylabel(ax, strrep(YVar,'_',' '), 'FontSize', options.AxFontSize);
        end
        title(ax, strrep(pDef.label,'_',' '), ...
              'FontSize', options.AxFontSize, 'FontWeight','normal', ...
              'Interpreter','none');

        % Store handles
        safeField = matlab.lang.makeValidName(pDef.label);
        axs.(safeField)   = ax;
        if ~isempty(st)
            stats.(safeField) = st;
        end
    end
end

end % ── end main function ──────────────────────────────────────────────────


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ── buildPanelDefs ────────────────────────────────────────────────────────
function panels = buildPanelDefs(c1Vals, allLevels, nRows, nCols)
%BUILDPANELDEFS  Build nRows x nCols cell array of panel definition structs.
%
% Each panel struct has:
%   .label      : string label for title / field name
%   .conditions : Kx2 cell array of {c1val, c2val} pairs to overlay
%   .isCombined : logical — true for last-column panels and all of last row
%                 Used so binned plots know to pool expCat rather than split.

nC1  = numel(c1Vals);
nLev = numel(allLevels);
panels = cell(nRows, nCols);

% Rows 1..nC1 : individual Cond1 values
for ri = 1:nC1
    c1 = c1Vals(ri);
    % Columns 1..nLev : single level  (NOT combined)
    for ci = 1:nLev
        lv = allLevels(ci);
        panels{ri,ci} = struct( ...
            'label',       sprintf('%s_%s', c1, lv), ...
            'conditions',  {{char(c1), char(lv)}}, ...
            'isCombined',  false );
    end
    % Last column: combined across all levels for this Cond1  (combined)
    combConds = cell(nLev, 2);
    for ci = 1:nLev
        combConds{ci,1} = char(c1);
        combConds{ci,2} = char(allLevels(ci));
    end
    levelStr = strjoin(allLevels, '+');
    panels{ri, nCols} = struct( ...
        'label',       sprintf('%s_%s', c1, levelStr), ...
        'conditions',  {combConds}, ...
        'isCombined',  true );
end

% Last row : cross-condition overlays  (all combined)
lastRow = nRows;
for ci = 1:nLev
    lv = allLevels(ci);
    crossConds = cell(nC1, 2);
    for ri = 1:nC1
        crossConds{ri,1} = char(c1Vals(ri));
        crossConds{ri,2} = char(lv);
    end
    c1Str = strjoin(c1Vals, '+');
    panels{lastRow, ci} = struct( ...
        'label',       sprintf('%s_%s', c1Str, lv), ...
        'conditions',  {crossConds}, ...
        'isCombined',  true );
end
% Last panel: all conditions  (combined)
allConds = cell(nC1*nLev, 2);
idx = 1;
for ri = 1:nC1
    for ci = 1:nLev
        allConds{idx,1} = char(c1Vals(ri));
        allConds{idx,2} = char(allLevels(ci));
        idx = idx + 1;
    end
end
panels{lastRow, nCols} = struct( ...
    'label',       'all', ...
    'conditions',  {allConds}, ...
    'isCombined',  true );
end


% ── resolveColorMap ───────────────────────────────────────────────────────
function cmap = resolveColorMap(c1Vals, allLevels, userColors)

defaultPalettes = {'blues','oranges','greens','purples', ...
                   'indigos','magentas','teals','redoranges'};
cmap = containers.Map('KeyType','char','ValueType','char');

paletteIdx = 0;
for i = 1:numel(c1Vals)
    c1 = char(c1Vals(i));
    for j = 1:numel(allLevels)
        c2  = char(allLevels(j));
        key = [c1, '_', c2]; % use '__' if conditions already have '_' in them

        if ~isempty(userColors) && isa(userColors, 'containers.Map') ...
                && isKey(userColors, key)
            cmap(key) = userColors(key);
        else
            paletteIdx = paletteIdx + 1;
            cmap(key)  = defaultPalettes{mod(paletteIdx-1, ...
                             numel(defaultPalettes)) + 1};
        end
    end
end
end


% ── collectConditionData ──────────────────────────────────────────────────
function condData = collectConditionData(tableT, XVar, YVar, ...
                        conditions, cond1Col, cond2Col, expCol)
%COLLECTCONDITIONDATA  For each (c1,c2) condition pair, collect per-ROI and
% pooled data.
%
%  roiX / roiY  : per-ROI means  (one value per ROI)
%  poolX / poolY: all raw values (one value per object inside each ROI)
%                 weighted so every ROI contributes equally (weight = 1/N
%                 for each of the N objects in that ROI).
%
% Both representations are always collected so the caller can choose which
% to use via ScatterMode.

nCond = size(conditions, 1);
condData(nCond) = struct( ...
    'label',     '', ...
    'c1',        '', ...
    'c2',        '', ...
    'roiX',      [], ...
    'roiY',      [], ...
    'expLabels', {{}}, ...       % per-ROI exp label
    'poolX',     [], ...
    'poolY',     [], ...
    'poolW',     [], ...
    'poolExpLabels', {{}} );     % per-object exp label (for 'all' mode)

for k = 1:nCond
    c1 = conditions{k,1};
    c2 = conditions{k,2};

    mask = string(tableT.(cond1Col)) == c1 & ...
           string(tableT.(cond2Col)) == c2;
    sub  = tableT(mask, :);
    nROI = height(sub);

    roiX   = nan(nROI, 1);
    roiY   = nan(nROI, 1);
    expLab = strings(nROI, 1);
    poolX  = [];
    poolY  = [];
    poolW  = [];
    poolExpLab = {};

    for ri = 1:nROI
        xVec = sub.(XVar){ri};
        yVec = sub.(YVar){ri};
        eLab = string(sub.(expCol)(ri));

        % per-ROI average
        roiX(ri) = mean(xVec, 'omitnan');
        roiY(ri) = mean(yVec, 'omitnan');
        expLab(ri) = eLab;

        % pooled with equal-ROI weighting
        n = numel(xVec);
        poolX = [poolX; xVec(:)];                        %#ok<AGROW>
        poolY = [poolY; yVec(:)];                        %#ok<AGROW>
        poolW = [poolW; repmat(1/n, n, 1)];              %#ok<AGROW>
        poolExpLab = [poolExpLab; repmat({char(eLab)}, n, 1)]; %#ok<AGROW>
    end

    % normalise weights so they sum to nROI (preserves scale)
    if ~isempty(poolW) && sum(poolW) > 0
        poolW = poolW / sum(poolW) * nROI;
    end

    condData(k).label          = sprintf('%s_%s', c1, c2);
    condData(k).c1             = c1;
    condData(k).c2             = c2;
    condData(k).roiX           = roiX;
    condData(k).roiY           = roiY;
    condData(k).expLabels      = cellstr(expLab);
    condData(k).poolX          = poolX;
    condData(k).poolY          = poolY;
    condData(k).poolW          = poolW;
    condData(k).poolExpLabels  = poolExpLab;
end
end


% ── plotScatter ───────────────────────────────────────────────────────────
function ellSt = plotScatter(ax, condData, colorMap, options, isCombined)
%PLOTSCATTER  Scatter plot coloured by expCat.
%
% ScatterMode = 'roi' : one dot per ROI (per-ROI mean)
% ScatterMode = 'all' : one dot per object (raw Nx1 values)
%
% ShowEllipse = true  : draw a 5%-95% percentile ellipse per condition
%                       (fit to whichever data mode is active).
% EllipseStats = true : when ≥2 conditions are present, test whether
%                       their covariance ellipses differ (Hotelling T²).
%                       Result is returned in ellSt.

ellSt = [];

useAll     = strcmpi(options.ScatterMode, 'all');
useHeatmap = strcmpi(options.ScatterStyle, 'heatmap');

for k = 1:numel(condData)
    cd   = condData(k);
    key  = [cd.c1, '_', cd.c2];
    pal  = colorMap(key);

    if useAll

        % --- Collect all points ---
        xAll = cd.poolX;
        yAll = cd.poolY;
        
        valid = isfinite(xAll) & isfinite(yAll);
        xAll = xAll(valid);
        yAll = yAll(valid);

        if useHeatmap
        
            if numel(xAll) < 10
                continue
            end
            
            % --- Compute density via 2D histogram ---
            nBins = 100;   % resolution (adjust)
            
            [N, Xedges, Yedges] = histcounts2(xAll, yAll, nBins);
            
            % Assign each point to a bin
            binX = discretize(xAll, Xedges);
            binY = discretize(yAll, Yedges);
            
            validBins = ~isnan(binX) & ~isnan(binY);
            linIdx = sub2ind(size(N), binX(validBins), binY(validBins));
            
            density = zeros(size(xAll));
            density(validBins) = N(linIdx);
            
            % Optional: log-scale (HIGHLY recommended)
            density = log1p(density);
            
            % --- Normalize density ---
            density = density / max(density);
            
            % --- Get condition colormap ---
            key = [cd.c1, '_', cd.c2];
            pal = colorMap(key);
            cmap = getPlotColors(pal, 256);
            
            % --- Map density → colors ---
            colorIdx = max(1, round(density * (size(cmap,1)-1)) + 1);
            colors = cmap(colorIdx, :);
            
            % --- Plot ---
            scatter(ax, xAll, yAll, options.MarkerSize, colors, ...
                'filled', ...
                'MarkerFaceAlpha', 0.8, ...
                'MarkerEdgeAlpha', 0.8);
            
            colormap(ax, cmap);
            colorbar(ax);
    
        else    
            % Group by expCat using poolExpLabels
            expGroups = unique(cd.poolExpLabels, 'stable');
            nExp      = numel(expGroups);
            colors    = getPlotColors(pal, nExp);
            markers   = getPlotMarkers(nExp);
    
            for e = 1:nExp
                eg  = expGroups{e};
                idx = strcmp(cd.poolExpLabels, eg);
                xE  = cd.poolX(idx);
                yE  = cd.poolY(idx);
                col = colors(e,:);
                mrk = markers{e};
    
                if options.ShowDots
                    scatter(ax, xE, yE, options.MarkerSize, ...
                        'Marker', mrk, ...
                        'MarkerFaceColor', col, ...
                        'MarkerEdgeColor', col, ...
                        'MarkerFaceAlpha', options.Alpha, ...
                        'MarkerEdgeAlpha', options.Alpha, ...
                        'DisplayName', sprintf('%s | %s', cd.label, eg));
                end
    
                if options.ShowMeans
                    scatter(ax, mean(xE,'omitnan'), mean(yE,'omitnan'), ...
                        options.MeanMarkerSize, ...
                        'Marker', mrk, ...
                        'MarkerFaceColor', col, ...
                        'MarkerEdgeColor', col*0.6, ...
                        'LineWidth', 1.2, ...
                        'DisplayName', sprintf('%s | %s (mean)', cd.label, eg));
                end
            end
        end

        % Ellipse drawn over all pooled data for this condition
        if options.ShowEllipse
            valid = isfinite(cd.poolX) & isfinite(cd.poolY);
            xAll  = cd.poolX(valid);
            yAll  = cd.poolY(valid);
            condCol = getPlotColors(pal, 1);   % representative colour
            drawPercentileEllipse(ax, xAll, yAll, condCol, 0.05, 0.95);
        end

    else
        % ── ROI mode: one dot per ROI ──────────────────────────────────
        expGroups = unique(cd.expLabels, 'stable');
        nExp      = numel(expGroups);
        colors    = getPlotColors(pal, nExp);
        markers   = getPlotMarkers(nExp);

        for e = 1:nExp
            eg  = expGroups{e};
            idx = strcmp(cd.expLabels, eg);
            xE  = cd.roiX(idx);
            yE  = cd.roiY(idx);
            col = colors(e,:);
            mrk = markers{e};

            if options.ShowDots
                scatter(ax, xE, yE, options.MarkerSize, ...
                    'Marker', mrk, ...
                    'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', col, ...
                    'MarkerFaceAlpha', options.Alpha, ...
                    'MarkerEdgeAlpha', options.Alpha, ...
                    'DisplayName', sprintf('%s | %s', cd.label, eg));
            end

            if options.ShowMeans
                scatter(ax, mean(xE,'omitnan'), mean(yE,'omitnan'), ...
                    options.MeanMarkerSize, ...
                    'Marker', mrk, ...
                    'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', col*0.6, ...
                    'LineWidth', 1.2, ...
                    'DisplayName', sprintf('%s | %s (mean)', cd.label, eg));
            end
        end

        % Ellipse drawn over all ROI-mean data for this condition
        if options.ShowEllipse
            valid = isfinite(cd.roiX) & isfinite(cd.roiY);
            xAll  = cd.roiX(valid);
            yAll  = cd.roiY(valid);
            condCol = getPlotColors(pal, 1);
            drawPercentileEllipse(ax, xAll, yAll, condCol, 0.05, 0.95);
        end
    end
end

% ── Ellipse comparison statistics (only when isCombined and ≥2 conds) ────
if options.ShowEllipse && options.EllipseStats && isCombined && numel(condData) >= 2
    ellSt = compareEllipses(condData, useAll);
    addEllipseStatsAnnotation(ax, ellSt, options.AxFontSize);
end
end


% ── drawPercentileEllipse ─────────────────────────────────────────────────
function drawPercentileEllipse(ax, x, y, col, pLow, pHigh)
%DRAWPERCENTILEELLIPSE  Draw an ellipse spanning the pLow–pHigh percentile
% range of the 2D data cloud.
%
% Method: fit a Gaussian to (x,y), then scale its covariance ellipse so
% that the semi-axes correspond to the requested quantiles of the marginal
% distributions.  This gives an intuitive "central region" oval that is
% robust against the number of points.

if numel(x) < 5, return; end

mu  = [mean(x,'omitnan'), mean(y,'omitnan')];
C   = cov(x, y);              % 2x2 covariance

% Eigen-decomposition for ellipse orientation
[V, D] = eig(C);
[~, si] = sort(diag(D));
V = V(:, si);
D = D(si, si);

% Scale semi-axes to requested percentile range (chi² for 2D Gaussian)
% The "radius" of the ellipse in SD units such that the area covers the
% requested fraction:  r² = chi2inv(pHigh - pLow, 2)
pFrac = pHigh - pLow;         % e.g. 0.90 for 5%-95%
% Use chi2inv approximation: chi2inv(p,2) = -2*ln(1-p)
r = sqrt(-2 * log(max(1 - pFrac, eps)));

theta   = linspace(0, 2*pi, 200);
circPts = r * [cos(theta); sin(theta)];   % unit circle scaled by r

% Stretch by sqrt of eigenvalues, rotate by eigenvectors
ell = V * sqrt(D) * circPts;

xEll = mu(1) + ell(1,:);
yEll = mu(2) + ell(2,:);

plot(ax, xEll, yEll, '-', ...
     'Color',     [col, 0.85], ...
     'LineWidth', 1.8, ...
     'DisplayName', sprintf('ellipse (%.0f%%–%.0f%%)', pLow*100, pHigh*100));
end


% ── compareEllipses ───────────────────────────────────────────────────────
function ellSt = compareEllipses(condData, useAll)
%COMPAREELLIPSES  Test whether the 2D distributions of conditions differ.
%
% For each pair of conditions, Hotelling's T² test is applied to the
% (x,y) cloud.  This tests whether the bivariate means differ, which is
% the parametric analogue of asking whether the ellipses are centred at
% the same location.
%
% Additionally, Box's M test is used to compare covariance matrices
% (i.e. whether the *shapes* of the ellipses differ).
%
% Returns a struct array with one element per pair:
%   .cond1, .cond2   : condition labels
%   .hotelling_T2    : T² statistic
%   .hotelling_F     : equivalent F statistic
%   .hotelling_p     : p-value (centre shift)
%   .boxM_chi2       : Box's M chi² approximation
%   .boxM_p          : p-value (shape/size difference)
%   .n1, .n2         : sample sizes

nCond = numel(condData);
pairs = nchoosek(1:nCond, 2);
ellSt(size(pairs,1)) = struct( ...
    'cond1','','cond2','', ...
    'hotelling_T2',NaN,'hotelling_F',NaN,'hotelling_p',NaN, ...
    'boxM_chi2',NaN,'boxM_p',NaN,'n1',0,'n2',0);

for pi = 1:size(pairs,1)
    i1 = pairs(pi,1);
    i2 = pairs(pi,2);
    cd1 = condData(i1);
    cd2 = condData(i2);

    if useAll
        v1 = isfinite(cd1.poolX) & isfinite(cd1.poolY);
        X1 = [cd1.poolX(v1), cd1.poolY(v1)];
        v2 = isfinite(cd2.poolX) & isfinite(cd2.poolY);
        X2 = [cd2.poolX(v2), cd2.poolY(v2)];
    else
        v1 = isfinite(cd1.roiX) & isfinite(cd1.roiY);
        X1 = [cd1.roiX(v1), cd1.roiY(v1)];
        v2 = isfinite(cd2.roiX) & isfinite(cd2.roiY);
        X2 = [cd2.roiX(v2), cd2.roiY(v2)];
    end

    ellSt(pi).cond1 = cd1.label;
    ellSt(pi).cond2 = cd2.label;
    ellSt(pi).n1    = size(X1,1);
    ellSt(pi).n2    = size(X2,1);

    if size(X1,1) < 5 || size(X2,1) < 5, continue; end

    % ── Hotelling's T² (centre shift) ─────────────────────────────────
    n1 = size(X1,1); n2 = size(X2,1); p = 2;
    mu1 = mean(X1); mu2 = mean(X2);
    S1  = cov(X1);  S2  = cov(X2);
    Sp  = ((n1-1)*S1 + (n2-1)*S2) / (n1+n2-2);   % pooled covariance

    if rcond(Sp) < eps
        % Degenerate covariance — skip
        continue
    end

    dm   = (mu1 - mu2)';
    T2   = (n1*n2/(n1+n2)) * (dm' / Sp * dm);
    df1  = p;
    df2  = n1 + n2 - p - 1;
    Fstat = T2 * (n1+n2-p-1) / (p*(n1+n2-2));
    pHot  = 1 - fcdf(Fstat, df1, df2);

    ellSt(pi).hotelling_T2 = T2;
    ellSt(pi).hotelling_F  = Fstat;
    ellSt(pi).hotelling_p  = pHot;

    % ── Box's M test (shape/size of ellipse) ──────────────────────────
    % Uses the chi² approximation (conservative but widely applicable).
    lnDet1 = logdet_safe(S1);
    lnDet2 = logdet_safe(S2);
    lnDetP = logdet_safe(Sp);

    M = (n1+n2-2)*lnDetP - (n1-1)*lnDet1 - (n2-1)*lnDet2;

    % chi² approximation with p*(p+1)/2 * (k-1) degrees of freedom (k=2 groups)
    df_box = p*(p+1)/2;
    pBox   = 1 - chi2cdf(M, df_box);

    ellSt(pi).boxM_chi2 = M;
    ellSt(pi).boxM_p    = pBox;
end
end


% ── addEllipseStatsAnnotation ─────────────────────────────────────────────
function addEllipseStatsAnnotation(ax, ellSt, fontSize)
%ADDELLIPSESTATSANNOTATION  Add ellipse comparison stats to the axes.
%
% Shows one line per pair:
%   "CondA vs CondB: centre p=X.XXX, shape p=X.XXX"

lines = cell(numel(ellSt), 1);
for k = 1:numel(ellSt)
    e = ellSt(k);
    lines{k} = sprintf('%s vs %s:  centre p=%s, shape p=%s', ...
        strrep(e.cond1,'_',' '), strrep(e.cond2,'_',' '), ...
        formatPval(e.hotelling_p), formatPval(e.boxM_p));
end
txt = strjoin(lines, newline);
text(ax, 0.03, 0.03, txt, ...
     'Units','normalized', ...
     'VerticalAlignment','bottom', ...
     'HorizontalAlignment','left', ...
     'FontSize', fontSize-1, ...
     'Interpreter','none', ...
     'BackgroundColor',[1 1 1 0.6], ...
     'Margin', 2);
end


% ── plotBinned ────────────────────────────────────────────────────────────
function axOut = plotBinned(ax, condData, colorMap, options, globalXLim, isCombined)
%PLOTBINNED  Bin XVar and show mean ± SEM ribbon of YVar per bin.
%
% Data used: all raw values (poolX/poolY with per-ROI equal weighting).
% This means if there are 20 objects per ROI you see 20 data points
% contributing to each bin — but each ROI contributes equally regardless
% of how many objects it contains.
%
% Panel behaviour:
%   isCombined = false  (top rows, individual level columns):
%       One ribbon per expCat group.
%
%   isCombined = true   (last column of top rows, entire last row):
%       All expCat pooled into a single ribbon per condition.
%       Colour = last (darkest) colour of the condition's palette.

if nargin < 5 || isempty(globalXLim) || ~all(isfinite(globalXLim))
    allX = [];
    for k = 1:numel(condData)
        allX = [allX; condData(k).poolX]; %#ok<AGROW>
    end
    valid = allX(isfinite(allX));
    globalXLim = [min(valid), max(valid)];
end
edges   = linspace(globalXLim(1), globalXLim(2), options.NBins+1);
binCtrs = (edges(1:end-1) + edges(2:end)) / 2;

axOut = ax;

for k = 1:numel(condData)
    cd  = condData(k);
    key = [cd.c1,'_',cd.c2];
    pal = colorMap(key);

    if isCombined
        % ── Single ribbon: pool all expCat, use last/darkest colour ───
        nColours = 6;   % request several shades, pick the darkest
        allColors = getPlotColors(pal, nColours);
        col = allColors(end,:);   % last = darkest shade

        valid = isfinite(cd.poolX) & isfinite(cd.poolY) & isfinite(cd.poolW);
        [mn, se] = weightedBinStats( ...
            cd.poolX(valid), cd.poolY(valid), cd.poolW(valid), edges);

        validBins = isfinite(mn) & isfinite(se);
        if sum(validBins) < 2, continue; end

        xV  = binCtrs(validBins);
        mnV = mn(validBins);
        seV = se(validBins);

        fill(ax, [xV, fliplr(xV)], ...
             [mnV+seV, fliplr(mnV-seV)], col, ...
             'FaceAlpha', 0.25, 'EdgeColor','none', ...
             'DisplayName', [cd.label, ' ±SEM']);
        plot(ax, xV, mnV, '-o', 'Color', col, ...
             'MarkerFaceColor', col, 'MarkerSize', 4, ...
             'LineWidth', 1.8, 'DisplayName', cd.label);

    else
        % ── One ribbon per expCat ──────────────────────────────────────
        expGroups = unique(cd.poolExpLabels, 'stable');
        nExp      = numel(expGroups);
        colors    = getPlotColors(pal, nExp);

        for e = 1:nExp
            eg   = expGroups{e};
            eMsk = strcmp(cd.poolExpLabels, eg);
            col  = colors(e,:);

            xE = cd.poolX(eMsk);
            yE = cd.poolY(eMsk);
            wE = cd.poolW(eMsk);

            valid = isfinite(xE) & isfinite(yE) & isfinite(wE);
            if sum(valid) < 2, continue; end

            [mn, se] = weightedBinStats(xE(valid), yE(valid), wE(valid), edges);

            validBins = isfinite(mn) & isfinite(se);
            if sum(validBins) < 2, continue; end

            xV  = binCtrs(validBins);
            mnV = mn(validBins);
            seV = se(validBins);

            fill(ax, [xV, fliplr(xV)], ...
                 [mnV+seV, fliplr(mnV-seV)], col, ...
                 'FaceAlpha', 0.25, 'EdgeColor','none', ...
                 'DisplayName', sprintf('%s | %s ±SEM', cd.label, eg));
            plot(ax, xV, mnV, '-o', 'Color', col, ...
                 'MarkerFaceColor', col, 'MarkerSize', 4, ...
                 'LineWidth', 1.8, ...
                 'DisplayName', sprintf('%s | %s', cd.label, eg));
        end
    end
end
end


% ── computeStats ──────────────────────────────────────────────────────────
function st = computeStats(condData, statTest, plotType, scatterMode, globalXLim, nBins)
%COMPUTESTATS  Compute weighted correlation per condition (pooled objects).

st = struct();
for k = 1:numel(condData)
    cd  = condData(k);
    key = matlab.lang.makeValidName([cd.c1,'_',cd.c2]);

    % N = ...
    if strcmpi(plotType, 'scatter') && strcmpi(scatterMode, 'roi')
        % --- ROI mode ---
        x = cd.roiX;
        y = cd.roiY;
        w = ones(size(x));   % equal weight per ROI
    else
        % --- pooled / "all" mode (default) ---
        x = cd.poolX;
        y = cd.poolY;
        w = cd.poolW;
    end

    valid = ~isnan(x) & ~isnan(y);
    x = x(valid); y = y(valid); w = w(valid);

    % restrict to points in bins that survived filtering
    if strcmpi(plotType, 'binned') && ~isempty(x)
        edges = linspace(globalXLim(1), globalXLim(2), nBins+1);
        keep  = getValidBinMask(x, y, w, edges);
        x = x(keep); y = y(keep); w = w(keep);
    end

    if numel(x) < 3
        st.(key) = struct('rho',NaN,'pval',NaN,'n',numel(x), ...
                          'test', statTest);
        continue
    end

    switch lower(statTest)
        case 'spearman'
            rx = weightedRank(x, w);
            ry = weightedRank(y, w);
            [rho, pval] = weightedCorr(rx, ry, w);
        case 'pearson'
            [rho, pval] = weightedCorr(x, y, w);
        otherwise
            error('Plot_CompareVar: unknown StatTest "%s"', statTest);
    end

    st.(key) = struct('rho',rho,'pval',pval,'n', numel(x), ...
                      'test', statTest);
end
end

% ── computeBinnedSlope ─────────────────────────────────────────────────────
function st = computeBinnedSlope(condData, globalXLim, nBins)
%COMPUTEBINNEDSLOPE  Weighted linear regression of binned means vs bin
% centres, per condition. Tests whether the ribbon's slope differs from
% zero — i.e. whether Y changes with X at a rate distinguishable from
% noise, using the same bins/filtering as the plotted ribbon.
%
% Weights = 1/se^2 per bin (more precise bins count more).
%
% Returns struct with one field per condition (key = c1_c2), each with:
%   .slope     : weighted LS slope (Y-units per X-unit)
%   .slope_se  : standard error of the slope
%   .pval      : p-value testing slope ~= 0 -> how different the slope is from zero -> < 0.05 -> confident slope is not zero
%   .n         : number of bins used in the fit
%   .intercept : weighted LS intercept (for reference)

edges   = linspace(globalXLim(1), globalXLim(2), nBins+1);
binCtrs = (edges(1:end-1) + edges(2:end)) / 2;

st = struct();
for k = 1:numel(condData)
    cd  = condData(k);
    key = matlab.lang.makeValidName([cd.c1,'_',cd.c2]);

    valid = isfinite(cd.poolX) & isfinite(cd.poolY) & isfinite(cd.poolW);
    [mn, se] = weightedBinStats( ...
        cd.poolX(valid), cd.poolY(valid), cd.poolW(valid), edges);

    validBins = isfinite(mn) & isfinite(se) & se > 0;
    n = sum(validBins);

    if n < 3
        st.(key) = struct('slope',NaN,'slope_se',NaN,'pval',NaN, ...
                          'n',n,'intercept',NaN);
        continue
    end

    x = binCtrs(validBins)';
    y = mn(validBins)';
    w = 1 ./ (se(validBins)'.^2);   % precision weighting

    W    = sum(w);
    xbar = sum(w .* x) / W;
    ybar = sum(w .* y) / W;
    Sxx  = sum(w .* (x - xbar).^2);
    Sxy  = sum(w .* (x - xbar) .* (y - ybar));

    if Sxx < eps
        st.(key) = struct('slope',NaN,'slope_se',NaN,'pval',NaN, ...
                          'n',n,'intercept',NaN);
        continue
    end

    b = Sxy / Sxx;
    a = ybar - b * xbar;

    resid = y - (a + b*x);
    SSE   = sum(w .* resid.^2);
    df    = n - 2;
    sigma2 = SSE / max(df, 1);
    se_b   = sqrt(sigma2 / Sxx);

    t    = b / se_b;
    pval = 2 * (1 - tcdf(abs(t), max(df,1)));

    st.(key) = struct('slope',b,'slope_se',se_b,'pval',pval, ...
                      'n',n,'intercept',a);
end
end

% ── addStatsAnnotation ────────────────────────────────────────────────────
function addStatsAnnotation(ax, st, fontSize)
%ADDSTATSANNOTATION  Add rho/p annotation to axes (top-left corner).

fields = fieldnames(st);
% Only print correlation sub-fields (skip ellipseTest struct)
lines  = {};
for k = 1:numel(fields)
    if strcmp(fields{k}, 'ellipseTest'), continue; end
    s = st.(fields{k});
    if ~isstruct(s) || ~isfield(s,'rho'), continue; end
    if isnan(s.rho)
        lines{end+1} = sprintf('%s: n.a.', fields{k}); %#ok<AGROW>
    else
        lines{end+1} = sprintf('\\rho=%.2f, p=%s (n=%d)', ...
            s.rho, formatPval(s.pval), s.n);            %#ok<AGROW>
    end
end
if isempty(lines), return; end
txt = strjoin(lines, newline);
text(ax, 0.03, 0.97, txt, ...
     'Units','normalized', ...
     'VerticalAlignment','top', ...
     'HorizontalAlignment','left', ...
     'FontSize', fontSize-1, ...
     'Interpreter','tex', ...
     'BackgroundColor',[1 1 1 0.6], ...
     'Margin', 2);
end

% ── addSlopeStatsAnnotation ────────────────────────────────────────────────
function addSlopeStatsAnnotation(ax, st, fontSize)
%ADDSLOPESTATSANNOTATION  Add binned-slope stats to the axes (top-right).

fields = fieldnames(st);
lines  = {};
for k = 1:numel(fields)
    s = st.(fields{k});
    if isnan(s.slope)
        lines{end+1} = sprintf('%s: slope n.a.', fields{k}); %#ok<AGROW>
    else
        lines{end+1} = sprintf('slope=%.3g, p=%s (bins=%d)', ...
            s.slope, formatPval(s.pval), s.n);              %#ok<AGROW>
    end
end
if isempty(lines), return; end
txt = strjoin(lines, newline);
text(ax, 0.97, 0.03, txt, ...
     'Units','normalized', ...
     'VerticalAlignment','bottom', ...
     'HorizontalAlignment','right', ...
     'FontSize', fontSize-1, ...
     'Interpreter','tex', ...
     'BackgroundColor',[1 1 1 0.6], ...
     'Margin', 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WEIGHTED STATISTICS HELPERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function validMask = getValidBinMask(x, ~, w, edges, minFrac, minEff)
%GETVALIDBINMASK  True for points belonging to a bin that has at least
% minFrac of total weight AND at least minEff effective ROIs.

if nargin < 5, minFrac = 0.02; end  % Discard bins holding less than 3% of total weight
if nargin < 6, minEff  = 15;   end  % “I want at least X independent, equally contributing ROIs in this bin”

binIdx    = discretize(x, edges);
nBins     = numel(edges) - 1;
totalW    = sum(w);
validMask = false(size(x));

for b = 1:nBins
    bMask = binIdx == b;
    if sum(bMask) < 2, continue; end
    wB = w(bMask);
    if sum(wB)/totalW < minFrac, continue; end
    nEff = 1 / sum((wB/sum(wB)).^2);
    if nEff < minEff, continue; end
    validMask(bMask) = true;
end
end

% ── weightedBinStats ──────────────────────────────────────────────────────
function [mn, se] = weightedBinStats(x, y, w, edges)
%WEIGHTEDBINSTATS  Weighted mean and SEM of y per x-bin.

nBins = numel(edges) - 1;
mn = nan(1, nBins);
se = nan(1, nBins);

binIdx = discretize(x, edges);
validMask = getValidBinMask(x, y, w, edges);

for b = 1:nBins
    bMask = binIdx == b & validMask;
    if sum(bMask) < 2, continue; end

    yB = y(bMask);
    wB = w(bMask);

    % Winsorization step to be less sensitive to outliers
    p = 5;   % percentage (5% / 95%)
    lo = prctile(yB, p);
    hi = prctile(yB, 100 - p);
    yB = max(min(yB, hi), lo);

    wB = wB / sum(wB);

    % mean
    mn(b) = sum(wB .* yB);
    % SEM
    se(b) = sqrt(sum(wB .* (yB - mn(b)).^2) / max(1 - sum(wB.^2), eps));
end

end


% ── weightedRank ──────────────────────────────────────────────────────────
function r = weightedRank(x, w)
%WEIGHTEDRANK  Weighted rank: midpoint of cumulative weight interval.

[~, si] = sort(x);
ws = w(si);
cr = cumsum(ws) - ws/2;
r  = zeros(size(x));
r(si) = cr;
end


% ── weightedCorr ──────────────────────────────────────────────────────────
function [rho, pval] = weightedCorr(x, y, w)
%WEIGHTEDCORR  Weighted Pearson correlation and approximate p-value.

w  = w / sum(w);
mx = sum(w .* x);
my = sum(w .* y);
sx = sqrt(sum(w .* (x - mx).^2));
sy = sqrt(sum(w .* (y - my).^2));

if sx < eps || sy < eps
    rho = NaN; pval = NaN; return
end

rho  = sum(w .* (x - mx) .* (y - my)) / (sx * sy);
rho  = max(-1, min(1, rho));

nEff = 1 / sum(w.^2);
t    = rho * sqrt(max(nEff-2,1)) / sqrt(max(1-rho^2, eps));
pval = 2 * (1 - tcdf(abs(t), max(nEff-2, 1)));
end


% ── logdet_safe ───────────────────────────────────────────────────────────
function ld = logdet_safe(A)
%LOGDET_SAFE  Log-determinant via Cholesky; falls back to eig for PSD mats.

try
    L  = chol(A);
    ld = 2 * sum(log(diag(L)));
catch
    ev = eig(A);
    ev = ev(ev > eps);
    ld = sum(log(ev));
end
end


% ── formatPval ────────────────────────────────────────────────────────────
function s = formatPval(p)
%FORMATPVAL  Format p-value for display.

if isnan(p),   s = 'n.a.'; return; end
if p < 0.001,  s = '<0.001'; return; end
s = sprintf('%.3f', p);
end