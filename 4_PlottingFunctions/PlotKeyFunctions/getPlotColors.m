function colors = getPlotColors(paletteName, n)
%GETPLOTCOLORS  Return an n-by-3 RGB color matrix for scatter plots.
%
%   colors = getPlotColors(paletteName, n)
%
% Inputs
%   paletteName : 'distinct' (default), 'lines', 'parula', 'turbo', 'hsv',
%                 'jet', 'hot', 'cool', 'spring', 'summer', 'autumn',
%                 'winter', 'gray', 'bone', 'copper', 'pink',
%                 'colorcube', 'flag', 'prism',
%                 'blues'   (light → dark blue)
%                 'oranges' (light → dark orange)
%                 'greens'  (light → dark green)
%                 'purples' (light → dark purple)
%
% Output
%   colors      : n-by-3 RGB matrix in [0,1]

if nargin < 1 || isempty(paletteName), paletteName = 'distinct'; end
if nargin < 2 || isempty(n),           n = 1;                    end
n = max(n, 1);

% ── Palettes backed by a continuous MATLAB colormap ──────────────────────
continuousMaps = {"parula","turbo","hsv","jet","hot","cool", ...
                  "spring","summer","autumn","winter",       ...
                  "gray","bone","copper","pink"};

switch lower(string(paletteName))

    % ── Categorical maps ──────────────────────────────────────────────────
    case "lines"
        colors = lines(n);

    case "colorcube"
        colors = colorcube(n);

    case "flag"
        base   = flag(max(n, 4));
        idx    = mod((0:n-1), size(base,1)) + 1;
        colors = base(idx, :);

    case "prism"
        base   = prism(max(n, 6));
        idx    = mod((0:n-1), size(base,1)) + 1;
        colors = base(idx, :);

    % ── Continuous maps: sample n evenly spaced rows ──────────────────────
    case continuousMaps
        base   = feval(lower(string(paletteName)), max(n, 256));
        idx    = round(linspace(1, size(base,1), n));
        colors = base(idx, :);

    % ── Single-hue shade palettes ─────────────────────────────────────────
    % RED
    case "redes"
        anchors = [0.404, 0.000, 0.051    % dark red / maroon
                   0.843, 0.188, 0.152;   % mid red
                   0.988, 0.573, 0.447;   % soft salmon
                   0.996, 0.878, 0.823];  % very light pink
        anchors = flipud(anchors); % go from light to dark
        colors = sampleAnchors(anchors, n);

    case "redmagentas"
        anchors = [0.980, 0.878, 0.878    % very light rose-red
                   0.937, 0.500, 0.550;   % soft crimson-pink
                   0.800, 0.150, 0.350;   % crimson-raspberry
                   0.550, 0.000, 0.200];  % deep red-magenta
        colors = sampleAnchors(anchors, n);

    case "magentas"
        anchors = [0.980, 0.921, 0.941    % very light blush
                   0.937, 0.627, 0.745;   % light rose
                   0.843, 0.188, 0.478;   % raspberry
                   0.490, 0.000, 0.235];  % deep magenta-red
        colors = sampleAnchors(anchors, n);

    case "vermillions"
        anchors = [0.996, 0.878, 0.843    % very light pinkish-orange
                   0.980, 0.580, 0.420;   % light vermillion
                   0.870, 0.220, 0.110;   % mid vermillion, red-dominant
                   0.550, 0.050, 0.020];  % dark, close to deep red
        colors = sampleAnchors(anchors, n);

    case "redoranges"
        anchors = [0.998, 0.929, 0.882    % very light peach
                   0.992, 0.733, 0.518;   % light coral
                   0.937, 0.396, 0.282;   % vermillion
                   0.647, 0.058, 0.082];  % dark red-orange
        colors = sampleAnchors(anchors, n);

    % BLUE
    case "blues"
        anchors = [0.031, 0.188, 0.420    % dark navy
                   0.129, 0.443, 0.710;   % mid blue
                   0.420, 0.682, 0.839;   % sky blue
                   0.870, 0.921, 0.969];  % very light blue
        anchors = flipud(anchors); % go from light to dark
        colors = sampleAnchors(anchors, n);

    case "blueindigos"
        anchors = [0.878, 0.914, 0.969    % very light blue-lavender
                   0.490, 0.620, 0.820;   % soft slate-blue
                   0.250, 0.380, 0.720;   % mid blue-indigo
                   0.080, 0.100, 0.500];  % deep blue-indigo
        colors = sampleAnchors(anchors, n);

    case "indigos"
        anchors = [0.925, 0.933, 0.980    % very light indigo
                   0.651, 0.741, 0.858;   % light steel blue
                   0.349, 0.467, 0.713;   % mid indigo
                   0.129, 0.188, 0.490];  % deep indigo
        colors = sampleAnchors(anchors, n);

    case "teals"
        anchors = [0.878, 0.949, 0.949    % very light blue-teal
                   0.502, 0.808, 0.820;   % light cyan-blue
                   0.098, 0.561, 0.620;   % mid teal, still blue-dominant
                   0.000, 0.278, 0.392];  % dark teal, close to navy
        colors = sampleAnchors(anchors, n);

    case "bluegreens"
        anchors = [0.878, 0.952, 0.941    % very light aqua
                   0.600, 0.847, 0.788;   % light teal
                   0.200, 0.627, 0.573;   % mid teal
                   0.000, 0.353, 0.392];  % dark teal
        colors = sampleAnchors(anchors, n);

    % OTHER
    case "yellows"
        anchors = [0.996, 0.973, 0.851    % very light cream-yellow
                   0.988, 0.902, 0.502;   % light yellow
                   0.961, 0.796, 0.051;   % your mid yellow
                   0.600, 0.471, 0.000];  % dark golden-yellow
        colors = sampleAnchors(anchors, n);

    case "oranges"
        anchors = [0.549, 0.180, 0.016    % dark burnt orange
                   0.890, 0.412, 0.078;   % mid orange
                   0.992, 0.682, 0.380;   % peach-orange
                   0.996, 0.902, 0.808];  % very light orange  
        anchors = flipud(anchors); % go from light to dark
        colors = sampleAnchors(anchors, n);

    case "greens"
        anchors = [0.000, 0.267, 0.106    % dark forest green
                   0.137, 0.545, 0.271;   % mid green
                   0.525, 0.808, 0.596;   % mint green
                   0.851, 0.941, 0.827];  % very light green  
        anchors = flipud(anchors); % go from light to dark
        colors = sampleAnchors(anchors, n);

    case "purples"
        anchors = [0.247, 0.000, 0.490    % deep violet
                   0.502, 0.212, 0.627;   % mid purple
                   0.737, 0.506, 0.741;   % mid lavender
                   0.937, 0.855, 0.961];  % very light purple/lavender  
        anchors = flipud(anchors); % go from light to dark
        colors = sampleAnchors(anchors, n);

    case "bluepurples"
        anchors = [0.878, 0.890, 0.973    % very light blue-lavender
                   0.557, 0.620, 0.949;   % light blue-purple
                   0.278, 0.376, 0.914;   % your mid blue-purple
                   0.098, 0.098, 0.549];  % dark blue-purple / deep indigo
        colors = sampleAnchors(anchors, n);

    % ── Single-hue shade palettes - classification────────────────────────
    case "bluepurple_custom"
        anchors = [ ...
            0.900, 0.910, 0.985    % very light cool lavender-blue
            0.520, 0.560, 0.900    % muted blue-purple
            0.220, 0.260, 0.780    % *** richer indigo (updated base feel)
            0.050, 0.060, 0.350];  % deep navy-purple
        colors = sampleAnchors(anchors, n);

    case "green_custom"
        anchors = [ ...
            0.851, 0.941, 0.827    % very light green
            0.525, 0.808, 0.596    % light green
            0.302, 0.686, 0.290    % YOUR base color
            0.000, 0.267, 0.106];  % dark forest green
        colors = sampleAnchors(anchors, n);

    case "yellow_custom"
            anchors = [ ...
                0.999, 0.985, 0.900    % very light lemon
                0.995, 0.940, 0.400    % bright light yellow
                0.980, 0.870, 0.000    % *** more pure yellow base
                0.650, 0.550, 0.000];  % darker mustard (less orange than before)
            colors = sampleAnchors(anchors, n);

    case "purple_custom"
        anchors = [ ...
            0.937, 0.855, 0.961    % very light lavender
            0.737, 0.506, 0.741    % light purple
            0.506, 0.278, 0.631    % YOUR base color
            0.247, 0.000, 0.490];  % deep violet
        colors = sampleAnchors(anchors, n);

    % ── Default: distinct evenly-spaced hues ─────────────────────────────
    otherwise
        h      = linspace(0, 1, n+1);
        h      = h(1:end-1);
        s      = 0.85;
        v      = 0.90;
        colors = hsv2rgb([h(:), repmat(s,n,1), repmat(v,n,1)]);
end
end


%% ── helper ───────────────────────────────────────────────────────────────
function colors = sampleAnchors(anchors, n)
%SAMPLEANCHORS  Linearly interpolate n colours from an anchor matrix.
    na = size(anchors, 1);

    if n == 1
        colors = anchors(round(0.25 * na), :); % avoid very light end
        return
    end

    % Avoid very light tail (e.g. use 15%–100% of the palette)
    tMin = 0.15;   % increase if still too light (e.g. 0.2–0.3)
    tMax = 1.0;

    t  = linspace(tMin, tMax, n);
    ta = linspace(0, 1, na);

    colors = zeros(n, 3);
    for c = 1:3
        colors(:,c) = interp1(ta, anchors(:,c), t, 'linear');
    end
end

% function colors = getPlotColors(paletteName, n)
% %GETPLOTCOLORS  Return an n-by-3 RGB color matrix for scatter plots.
% %
% %   colors = getPlotColors(paletteName, n)
% %
% % Inputs
% %   paletteName : 'distinct' (default), 'lines', 'parula', 'turbo', 'hsv'
% %   n           : number of colors needed
% %
% % Output
% %   colors      : n-by-3 RGB matrix in [0,1]
% 
% if nargin < 1 || isempty(paletteName), paletteName = 'distinct'; end
% if nargin < 2 || isempty(n),           n = 1;                    end
% 
% n = max(n, 1);
% 
% switch lower(string(paletteName))
%     case "lines"
%         colors = lines(n);
% 
%     case "parula"
%         big    = parula(max(n, 256));
%         idx    = round(linspace(1, size(big,1), n));
%         colors = big(idx, :);
% 
%     case "turbo"
%         big    = turbo(max(n, 256));
%         idx    = round(linspace(1, size(big,1), n));
%         colors = big(idx, :);
% 
%     case "hsv"
%         colors = hsv(n);
% 
%     otherwise  % 'distinct': evenly spaced hues, high saturation/value
%         h      = linspace(0, 1, n+1);
%         h      = h(1:end-1);           % drop duplicate at 1
%         s      = 0.85;
%         v      = 0.90;
%         colors = hsv2rgb([h(:), repmat(s, n, 1), repmat(v, n, 1)]);
% end
% end
