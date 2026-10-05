function fig = plotScatter2var(T, VarX, VarY, varargin)
% plotScatter2var
% Creates a 3x5 grid of scatter plots:
% Row 1: CT only, Row 2: CCM only, Row 3: CT vs CCM combined per extrainfoCat.
% Colors & markers in rows 1-2 = DataFolder; Row 3 = CT (blue) vs CCM (red) + regression lines.
%
% REQUIRED columns in T:
%   VarX, VarY, disease (CT/CCM), extrainfoCat, DataFolder
%
% Name-Value options:
%   'Columns'    : extrainfoCat labels for 5 columns (default: {"rest","R1","R2","10kPa","1kPa"})
%   'Palette'    : "distinct"|"lines"|"parula"|"hsv"|"turbo" (default: "distinct")
%   'MarkerSize' : scalar (default: 36)
%   'Alpha'      : 0..1 (default: 0.85)
%   'XCutoff'    : numeric (default: -Inf)
%   'YCutoff'    : numeric (default: -Inf)
%   'Title'      : main title
%
% Example:
% fig = plotScatter2var(T, 'XCutoff',100, 'Palette',"turbo");

% ----------------------- Parse inputs -----------------------
p = inputParser;
p.addParameter('includePatterns', {}, @(c) iscellstr(c) || isstring(c));
p.addParameter('excludePatterns', {}, @(c) iscellstr(c) || isstring(c));
p.addParameter('Columns', {"rest","R1","R2","10kPa","1kPa"}, @(c) iscellstr(c) || isstring(c));
p.addParameter('Palette', "distinct", @(s) isstring(s) || ischar(s));
p.addParameter('MarkerSize', 36, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('Alpha', 0.85, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('XCutoff', -Inf, @(x) isnumeric(x) && isscalar(x));
p.addParameter('YCutoff', -Inf, @(x) isnumeric(x) && isscalar(x));
p.addParameter('Title', sprintf('%s vs %s — CT/CCM by extrainfoCat', VarX, VarY), @(s) ischar(s) || isstring(s));
p.parse(varargin{:});

includePatterns = string(p.Results.includePatterns);
excludePatterns = string(p.Results.excludePatterns);
cols        = string(p.Results.Columns);
paletteName = string(p.Results.Palette);
mSize       = p.Results.MarkerSize;
alphaVal    = p.Results.Alpha;
xCutoff     = p.Results.XCutoff;
yCutoff     = p.Results.YCutoff;
mainTitle   = char(p.Results.Title);

% ----------------------- Filter columns -------------------
T = Filter_IncludeExclude(T, includePatterns, excludePatterns);

% ----------------------- Validate columns -------------------
requiredVars = {VarX, VarY,'disease','extrainfoCat','DataFolder'};
missingBase  = setdiff(requiredVars, T.Properties.VariableNames);
if ~isempty(missingBase)
    error('Missing required column(s): %s', strjoin(missingBase, ', '));
end

x = double(T.(VarX));
y = double(T.(VarY));
diseaseStr = upper(string(T.disease));
extra = string(T.extrainfoCat);
%extra = categorical(string(T.extrainfoCat));%, cols, 'Ordinal', true);
folder = string(T.DataFolder);

validRows = ~isnan(x) & ~isnan(y) & (diseaseStr=="CT" | diseaseStr=="CCM") ...
            & ~ismissing(extra) & ~ismissing(folder) ...
            & (x >= xCutoff) & (y >= yCutoff);
if ~any(validRows)
    warning('No valid rows after filtering.');
    fig = [];
    return;
end

xF = x(validRows); yF = y(validRows);
diseaseF = diseaseStr(validRows);
extraF = extra(validRows);
folderF = folder(validRows);

% ----------------------- Colors & markers -------------------
folderCats = unique(folderF, 'stable');
nF = numel(folderCats);
colorsF = selectPaletteColors(paletteName, nF);
markerSet = {'o','s','^','d','v','>','<','p','h','x','+','*'};
markersF = markerSet(mod(0:nF-1, numel(markerSet)) + 1);

% Axis limits
xLimAll = [min(xF)-0.05*range(xF), max(xF)+0.05*range(xF)];
yLimAll = [min(yF)-0.05*range(yF), max(yF)+0.05*range(yF)];

% ----------------------- Figure & layout --------------------
fig = figure('Color','w');
tl = tiledlayout(fig,3,6,'TileSpacing','compact','Padding','compact');
try, title(tl, mainTitle, 'Interpreter','none'); catch, sgtitle(mainTitle); end

lgHandles = gobjects(nF,1); hasHandle = false(nF,1);
rowNames = ["CT","CCM"];

% ---------------- Row 1 (CT) ----------------
for c = 1:5
    ax = nexttile(tl, c);                    % tiles 1..5
    hold(ax,'on'); grid(ax,'on'); box(ax,'off');
    thisCat = cols(c);
    title(ax, sprintf('CT — %s', thisCat), 'Interpreter','none');

    idxBase = (diseaseF=="CT") & strcmpi(extraF, thisCat);
    if any(idxBase)
        for f = 1:nF
            idx = idxBase & strcmp(folderF,folderCats(f));
            if ~any(idx), continue; end
            h = scatter(ax,xF(idx),yF(idx),mSize,...
                'Marker',markersF{f},...
                'MarkerFaceColor',colorsF(f,:),...
                'MarkerEdgeColor',colorsF(f,:),...
                'MarkerFaceAlpha',alphaVal,'MarkerEdgeAlpha',alphaVal,...
                'DisplayName',char(folderCats(f)));
            if ~hasHandle(f), lgHandles(f)=h; hasHandle(f)=true; end
        end
    end

    xlabel(ax, VarX, 'Interpreter', 'none'); ylabel(ax, VarY, 'Interpreter', 'none');
    xlim(ax, xLimAll); ylim(ax, yLimAll);
    %if isfinite(xCutoff), xline(ax, xCutoff, 'k:'); end
    %if isfinite(yCutoff), yline(ax, yCutoff, 'k:'); end
end

% ---------------- Row 2 (CCM) ----------------
for c = 1:5
    ax = nexttile(tl, 6 + c);               % tiles 7..11
    hold(ax,'on'); grid(ax,'on'); box(ax,'off');
    thisCat = cols(c);
    title(ax, sprintf('CCM — %s', thisCat), 'Interpreter','none');

    idxBase = (diseaseF=="CCM") & strcmpi(extraF, thisCat);
    if any(idxBase)
        for f = 1:nF
            idx = idxBase & strcmp(folderF,folderCats(f));
            if ~any(idx), continue; end
            h = scatter(ax,xF(idx),yF(idx),mSize,...
                'Marker',markersF{f},...
                'MarkerFaceColor',colorsF(f,:),...
                'MarkerEdgeColor',colorsF(f,:),...
                'MarkerFaceAlpha',alphaVal,'MarkerEdgeAlpha',alphaVal,...
                'DisplayName',char(folderCats(f)));
            if ~hasHandle(f), lgHandles(f)=h; hasHandle(f)=true; end
        end
    end

    xlabel(ax, VarX, 'Interpreter', 'none'); ylabel(ax, VarY,'Interpreter', 'none');
    xlim(ax, xLimAll); ylim(ax, yLimAll);
    %if isfinite(xCutoff), xline(ax, xCutoff, 'k:'); end
    %if isfinite(yCutoff), yline(ax, yCutoff, 'k:'); end
end

% ----------------------- Row 3 (combined CT vs CCM) --------
colCT=[0 0.4470 0.7410]; colCCM=[0.8500 0.3250 0.0980];

for c = 1:5
    ax = nexttile(tl, 12 + c);   % <-- row 3, columns 1..5 => tiles 13..17
    hold(ax,'on'); grid(ax,'on'); box(ax,'off');

    thisCat = cols(c);
    title(ax, sprintf('CT vs CCM — %s', thisCat), 'Interpreter','none');

    idxCT  = (diseaseF=="CT")  & strcmpi(extraF, thisCat);
    idxCCM = (diseaseF=="CCM") & strcmpi(extraF, thisCat);

    if any(idxCT)
        scatter(ax, xF(idxCT),  yF(idxCT), mSize, 'o', ...
            'MarkerFaceColor', colCT, 'MarkerEdgeColor', colCT, ...
            'MarkerFaceAlpha', alphaVal);
        pCT = polyfit(xF(idxCT), yF(idxCT), 1);
        xx  = linspace(min(xF(idxCT)), max(xF(idxCT)), 100);
        plot(ax, xx, polyval(pCT, xx), 'Color', colCT, 'LineWidth', 1.5);
    end

    if any(idxCCM)
        scatter(ax, xF(idxCCM), yF(idxCCM), mSize, 'o', ...
            'MarkerFaceColor', colCCM, 'MarkerEdgeColor', colCCM, ...
            'MarkerFaceAlpha', alphaVal);
        pCCM = polyfit(xF(idxCCM), yF(idxCCM), 1);
        xx   = linspace(min(xF(idxCCM)), max(xF(idxCCM)), 100);
        plot(ax, xx, polyval(pCCM, xx), 'Color', colCCM, 'LineWidth', 1.5);
    end

    xlabel(ax, VarX, 'Interpreter', 'none'); ylabel(ax, VarY, 'Interpreter', 'none');
    xlim(ax, xLimAll); ylim(ax, yLimAll);
    %if isfinite(xCutoff), xline(ax, xCutoff, 'k:'); end
    %if isfinite(yCutoff), yline(ax, yCutoff, 'k:'); end
end

% ----------------------- Legend outside ---------------------
% Legend in row 1, col 6 (tile 6)
axLegend = nexttile(tl, 6);
cla(axLegend); axis(axLegend, 'off');  % keep blank
useIdx = find(hasHandle);
if ~isempty(useIdx)
    % Attach legend to the legend tile explicitly
    lg = legend(axLegend, lgHandles(useIdx), cellstr(folderCats(useIdx)), ...
        'Location', 'northwest', 'Interpreter', 'none');
    %try, title(lg, 'DataFolder'); end
end

% Force blanks in (row 2, col 6) and (row 3, col 6)
axBlank2 = nexttile(tl, 12);  % (2,6)
cla(axBlank2); axis(axBlank2, 'off');

axBlank3 = nexttile(tl, 18);  % (3,6)
cla(axBlank3); axis(axBlank3, 'off');

end

% ======================= Helpers ==============================
function colors=selectPaletteColors(paletteName,nF)
switch lower(paletteName)
    case "lines", colors=lines(nF);
    case "parula", big=parula(max(nF,256)); idx=round(linspace(1,size(big,1),nF)); colors=big(idx,:);
    case "turbo", big=turbo(max(nF,256)); idx=round(linspace(1,size(big,1),nF)); colors=big(idx,:);
    case "hsv", colors=hsv(nF);
    otherwise, colors=makeDistinctColors(nF);
end
end

function C=makeDistinctColors(n)
h=linspace(0,1,n+1); h=h(1:end-1);
s=0.85; v=0.9;
C=hsv2rgb([h(:),repmat(s,n,1),repmat(v,n,1)]);
end
