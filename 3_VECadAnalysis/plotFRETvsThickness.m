function fig = plotFRETvsThickness(T, varargin)
% plotFRETvsThickness
% Bins skel_thickness and plots mean FRET per bin,
% organized in a 3x6 grid matching plotScatter2var layout.
%
% REQUIRED columns in T:
%   skel_thickness, skel_FRET, disease (CT/CCM), extrainfoCat, DataFolder
%
% Name-Value options:
%   'Columns'    : extrainfoCat labels (default: {"rest","R1","R2","10kPa","1kPa"})
%   'NBins'      : number of thickness bins (default: 10)
%   'Palette'    : "distinct"|"lines"|"parula"|"hsv"|"turbo" (default: "distinct")
%   'Alpha'      : 0..1 (default: 0.85)
%   'XCutoff'    : max thickness to include in µm (default: Inf)
%   'Title'      : main figure title

% ----------------------- Parse inputs -----------------------
p = inputParser;
p.addParameter('XVar', 'thickness', @(s) ischar(s)||isstring(s));
p.addParameter('YVar', 'fret',      @(s) ischar(s)||isstring(s));
p.addParameter('includePatterns', {}, @(c) iscellstr(c) || isstring(c));
p.addParameter('excludePatterns', {}, @(c) iscellstr(c) || isstring(c));
p.addParameter('Columns',  {"rest","R1","R2","10kPa","1kPa"}, @(c) iscellstr(c)||isstring(c));
p.addParameter('NBins',    10,         @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('Palette',  "distinct", @(s) isstring(s)||ischar(s));
p.addParameter('Alpha',    0.85,       @(x) isnumeric(x) && x>=0 && x<=1);
p.addParameter('XCutoff',  Inf,        @(x) isnumeric(x) && isscalar(x));
p.addParameter('Title',    'Mean FRET vs Junction Thickness', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});

xVar            = lower(string(p.Results.XVar));  % 'thickness' | 'fret' | 'int'
yVar            = lower(string(p.Results.YVar));  % 'thickness' | 'fret' | 'int'
cols            = string(p.Results.Columns);
nBins           = p.Results.NBins;
paletteName     = string(p.Results.Palette);
alphaVal        = p.Results.Alpha;
xCutoff         = p.Results.XCutoff;
mainTitle       = char(p.Results.Title);
includePatterns = string(p.Results.includePatterns);
excludePatterns = string(p.Results.excludePatterns);

% ----------------------- Filter -----------------------
%T = Filter_IncludeExclude(T, includePatterns, excludePatterns);
if ~isempty(includePatterns) || ~isempty(excludePatterns)
    T = Filter_IncludeExclude(T, includePatterns, excludePatterns);
end

% ----------------------- Validate columns -----------------------
requiredVars = {'Thickness_Skeleton','Skel_FRET','Skel_Int','disease','extrainfoCat','DataFolder'};
missingVars  = setdiff(requiredVars, T.Properties.VariableNames);
if ~isempty(missingVars)
    error('Missing required column(s): %s', strjoin(missingVars, ', '));
end

% ----------------------- Fast unpack cell arrays -----------------------
% Remove empty rows first
hasData = ~cellfun(@isempty, T.Thickness_Skeleton);
T       = T(hasData, :);

nPerRow      = cellfun(@numel, T.Thickness_Skeleton);
allThickness = vertcat(T.Thickness_Skeleton{:});
allFRET      = vertcat(T.Skel_FRET{:});
allInt       = vertcat(T.Skel_Int{:});
allDisease   = repelem(upper(string(T.disease)),   nPerRow(:));
allExtra     = repelem(string(T.extrainfoCat),     nPerRow(:));
allFolder    = repelem(string(T.DataFolder),       nPerRow(:));

% Helper to to resolve variable name → data array
function vec = resolveVar(name, allThickness, allFRET, allInt)
    switch lower(string(name))
        case 'thickness', vec = allThickness;
        case 'fret',      vec = allFRET;
        case 'int',       vec = allInt;
        otherwise, error('Unknown variable "%s". Use thickness/fret/int.', name);
    end
end

xData = resolveVar(xVar, allThickness, allFRET, allInt);
yData = resolveVar(yVar, allThickness, allFRET, allInt);

% ----------------------- Global valid mask & bin edges -----------------------
validAll = xData > 0 & xData <= xCutoff & ~isnan(yData) & ~isnan(xData);

binEdges = linspace(min(xData(validAll)), ...
    min(max(xData(validAll)), xCutoff), nBins+1);
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

% Compute bin IDs for ALL data once
[~, ~, binIDAll] = histcounts(xData, binEdges);

% ----------------------- Colors per DataFolder -----------------------
folderCats = unique(allFolder(validAll), 'stable');
nF         = numel(folderCats);
colorsF    = selectPaletteColors(paletteName, nF);

% ----------------------- Pre-compute global y-limits -----------------------
yAll = [];
for c = 1:numel(cols)
    for f = 1:nF
        for d = ["CT","CCM"]
            idx = validAll & strcmpi(allExtra, cols(c)) & ...
                  allDisease == d & allFolder == folderCats(f);
            if ~any(idx), continue; end
            [bm, ~] = fastBinStats(yData(idx), binIDAll(idx), nBins);
            yAll    = [yAll; bm(~isnan(bm))];
        end
    end
end

if isempty(yAll), warning('No valid data to plot.'); fig = []; return; end
yLimAll = [min(yAll) - 0.05*range(yAll), max(yAll) + 0.05*range(yAll)];
xLimAll = [binEdges(1), binEdges(end)];

% testpatch
yPad = 0.05 * range(yAll);
if yPad == 0, yPad = 0.5; end          % fallback padding if all values identical
yLimAll = [min(yAll) - yPad, max(yAll) + yPad];
xPad = 0.05 * (binEdges(end) - binEdges(1));
if xPad == 0, xPad = 0.5; end
xLimAll = [binEdges(1) - xPad, binEdges(end) + xPad];

% ----------------------- Figure & layout -----------------------
fig = figure('Color','w');
tl  = tiledlayout(fig, 3, 6, 'TileSpacing','compact', 'Padding','compact');
try, title(tl, mainTitle, 'Interpreter','none'); catch, sgtitle(mainTitle); end

lgHandles = gobjects(nF, 1);
hasHandle = false(nF, 1);

% ----------------------- Rows 1 (CT) and 2 (CCM) -----------------------
for row = 1:2
    dLabel = ["CT","CCM"];
    d      = dLabel(row);

    for c = 1:5
        ax = nexttile(tl, (row-1)*6 + c);
        hold(ax,'on'); grid(ax,'on'); box(ax,'off');
        title(ax, sprintf('%s — %s', d, cols(c)), 'Interpreter','none');

        for f = 1:nF
            idx = validAll & strcmpi(allExtra, cols(c)) & ...
                  allDisease == d & allFolder == folderCats(f);
            if ~any(idx), continue; end

            [binMean, binSEM] = fastBinStats(yData(idx), binIDAll(idx), nBins);

            hasData2 = ~isnan(binMean);
            if ~any(hasData2), continue; end

            xp  = binCenters(hasData2);
            yp  = binMean(hasData2);
            ye  = binSEM(hasData2);
            col = colorsF(f,:);

            fill(ax, [xp(:)', fliplr(xp(:)')], [yp(:)'+ye(:)', fliplr(yp(:)'-ye(:)')], col, ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none');
            h = plot(ax, xp, yp, '-o', 'Color', col, 'LineWidth', 1.5, ...
                'MarkerFaceColor', col, 'MarkerSize', 5, ...
                'DisplayName', char(folderCats(f)));

            if ~hasHandle(f), lgHandles(f) = h; hasHandle(f) = true; end
        end

        xlabel(ax, xVar, 'Interpreter','none');
        ylabel(ax, yVar, 'Interpreter','none');
        xlim(ax, xLimAll); ylim(ax, yLimAll);

    end
end

% ----------------------- Row 3: CT vs CCM overlaid -----------------------
colCT  = [0.0000 0.4470 0.7410];
colCCM = [0.8500 0.3250 0.0980];
dColors = [colCT; colCCM];

for c = 1:5
    ax = nexttile(tl, 12 + c);
    hold(ax,'on'); grid(ax,'on'); box(ax,'off');
    title(ax, sprintf('CT vs CCM — %s', cols(c)), 'Interpreter','none');

    for di = 1:2
        d   = ["CT","CCM"]; 
        idx = validAll & strcmpi(allExtra, cols(c)) & allDisease == d(di);
        if ~any(idx), continue; end

        [binMean, binSEM] = fastBinStats(yData(idx), binIDAll(idx), nBins);

        hasData2 = ~isnan(binMean);
        if ~any(hasData2), continue; end

        xp  = binCenters(hasData2);
        yp  = binMean(hasData2);
        ye  = binSEM(hasData2);
        col = dColors(di,:);

        fill(ax, [xp(:)', fliplr(xp(:)')], [yp(:)'+ye(:)', fliplr(yp(:)'-ye(:)')], col, ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(ax, xp, yp, '-o', 'Color', col, 'LineWidth', 1.5, ...
            'MarkerFaceColor', col, 'MarkerSize', 5, ...
            'DisplayName', char(d(di)));
    end

    legend(ax, 'Location','best', 'Interpreter','none');
    xlabel(ax, xVar, 'Interpreter','none');
    ylabel(ax, yVar, 'Interpreter','none');

    xlim(ax, xLimAll); ylim(ax, yLimAll);
end

% ----------------------- Legend tile (row 1, col 6) -----------------------
axLegend = nexttile(tl, 6); cla(axLegend); axis(axLegend,'off');
useIdx = find(hasHandle);
if ~isempty(useIdx)
    legend(axLegend, lgHandles(useIdx), cellstr(folderCats(useIdx)), ...
        'Location','northwest', 'Interpreter','none');
end
nexttile(tl, 12); cla; axis off;
nexttile(tl, 18); cla; axis off;

end

% =====================================================================
% Helper: vectorized bin mean + SEM using accumarray
% =====================================================================
function [binMean, binSEM] = fastBinStats(fretVals, binIDs, nBins)
    binMean = nan(nBins, 1);
    binSEM  = nan(nBins, 1);

    % Only use in-range bins (histcounts returns 0 for out-of-range)
    ok      = binIDs > 0 & binIDs <= nBins;
    if ~any(ok), return; end

    ids   = binIDs(ok);
    vals  = double(fretVals(ok));

    counts = accumarray(ids, 1,        [nBins 1]);
    sums   = accumarray(ids, vals,     [nBins 1]);
    sumsq  = accumarray(ids, vals.^2,  [nBins 1]);

    hasEnough = counts >= 2;
    binMean(hasEnough) = sums(hasEnough) ./ counts(hasEnough);
    variance           = (sumsq(hasEnough) - sums(hasEnough).^2 ./ counts(hasEnough)) ...
                         ./ (counts(hasEnough) - 1);
    binSEM(hasEnough)  = sqrt(variance) ./ sqrt(counts(hasEnough));
end

% =====================================================================
% Helper: palette selection
% =====================================================================
function colors = selectPaletteColors(paletteName, nF)
    switch lower(paletteName)
        case "lines",  colors = lines(nF);
        case "parula", big=parula(max(nF,256)); idx=round(linspace(1,size(big,1),nF)); colors=big(idx,:);
        case "turbo",  big=turbo(max(nF,256));  idx=round(linspace(1,size(big,1),nF)); colors=big(idx,:);
        case "hsv",    colors = hsv(nF);
        otherwise,     colors = makeDistinctColors(nF);
    end
end


function C = makeDistinctColors(n)
    h = linspace(0,1,n+1); h = h(1:end-1);
    C = hsv2rgb([h(:), repmat(0.85,n,1), repmat(0.9,n,1)]);
end