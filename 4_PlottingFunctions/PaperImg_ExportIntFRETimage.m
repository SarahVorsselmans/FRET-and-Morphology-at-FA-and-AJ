function PaperImg_ExportIntFRETimage(T, Input, SampleFolder, EntryNumber, options)
% ExportIntFRETimage  Export Int and/or FRET images as .tif for a given
%                    sample/entry combination found in table T.
%
% Usage:
%   ExportIntFRETimage(T, Input, SampleFolder, EntryNumber)
%   ExportIntFRETimage(T, Input, SampleFolder, EntryNumber, options)
%
% Inputs:
%   T             - Table containing at least columns: SampleFolder,
%                   EntryNumber, Int, FRET (cell arrays of 512x512 doubles)
%   Input         - Struct with field OutPath (output directory)
%   SampleFolder  - String/char, value to match in T.SampleFolder
%   EntryNumber   - Scalar, value to match in T.EntryNumber
%
% Optional name-value inputs (options):
%   options.Int   - Logical, export Int image  (default: true)
%   options.FRET  - Logical, export FRET image (default: true)
%
% Output:
%   Figures are displayed and saved to Input.OutPath as:
%       <SampleFolder>_Cell_<EntryNumber>_Int.tif
%       <SampleFolder>_Cell_<EntryNumber>_FRET.tif

    arguments
        T               table
        Input           struct
        SampleFolder    {mustBeTextScalar}
        EntryNumber     (1,1) double
        options.Int       (1,1) logical = true
        options.FRET      (1,1) logical = true
        options.Class     (1,1) logical = false
        options.Skel      (1,1) logical = false
        options.IntLow    (1,1) double  = 1          % percentile for black point
        options.IntHigh   (1,1) double  = 99         % percentile for white point
        options.ThickMax  (1,1) double  = inf  
        % options.ScaleBar     (1,1) logical = true    % toggle scale bar on/off
        options.ScaleBarUm   (1,1) double  = 20      % scale bar length in µm
        options.ExportSVG    (1,1) logical = false   % default off, since .tif is primary
    end

    %% --- Locate the unique row in T -----------------------------------
    rowMask = strcmp(T.SampleFolder, SampleFolder) & ...
              (T.EntryNumber == EntryNumber);

    nHits = sum(rowMask);
    if nHits == 0
        error('ExportIntFRETimage:noMatch', ...
              'No row found for SampleFolder="%s" and EntryNumber=%d.', ...
              SampleFolder, EntryNumber);
    elseif nHits > 1
        error('ExportIntFRETimage:ambiguousMatch', ...
              'Multiple rows (%d) found for SampleFolder="%s" and EntryNumber=%d.', ...
              nHits, SampleFolder, EntryNumber);
    end

    rowIdx   = find(rowMask);
    outPath  = Input.OutPath;
    baseName = sprintf('%s_Cell_%d', SampleFolder, EntryNumber);

    %% --- Intensity image ----------------------------------------------
    if options.Int
        A   = T.IntD{rowIdx};
        low  = prctile(A(:), options.IntLow);
        high = prctile(A(:), options.IntHigh);

        figInt = figure('Name', [baseName ' – Intensity']);
        % imagesc(A, [O 200]);
        imagesc(A, [low high]);
        colormap(figInt, gray);
        colorbar;
        axis square;
        axis off;
        title(strrep(baseName, '_', '\_'), 'FontWeight', 'normal');

        % --- Scale bar ---
        pixSize      = Input.PixSize;          % µm per pixel
        barLengthUm  = options.ScaleBarUm;     % µm
        barLengthPx  = barLengthUm / pixSize;  % pixels
        % Position: bottom-right corner with a small margin
        dimension = Input.Image_Resolution;
        margin    = 20;                        % pixels from edge
        barX      = dimension - margin - barLengthPx;
        barY      = dimension - margin;
        barThick  = 8;                         % bar height in pixels
        rectangle('Position', [barX, barY, barLengthPx, barThick], ...
                  'FaceColor', 'white', 'EdgeColor', 'white');
        text(barX + barLengthPx/2, barY - 5, '20 µm', ...
             'Color', 'white', 'FontSize', 10, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

        % Save
        savePath = fullfile(outPath, [baseName '_Int.tif']);
        exportgraphics(figInt, savePath, 'Resolution', 300);
        fprintf('Saved Int image: %s\n', savePath);
        if options.ExportSVG
            print(figInt, [fullfile(outPath, [baseName '_Int']) '.svg'], '-dsvg');
        end
    end

    %% --- FRET image ---------------------------------------------------
    if options.FRET
        BlackJet = loadBlackJet();   % helper nested function below

        F       = T.FRETimageFA{rowIdx};

        figFRET = figure('Name', [baseName ' – FRET']);
        imagesc(F, [0 40]);
        colormap(figFRET, BlackJet);
        c = colorbar;
        c.Label.String = 'Apparent FRET Efficiency %';
        axis square;
        axis off;
        title(strrep(baseName, '_', '\_'), 'FontWeight', 'normal');

        % save
        savePath = fullfile(outPath, [baseName '_FRET.tif']);
        exportgraphics(figFRET, savePath, 'Resolution', 300);
        fprintf('Saved FRET image: %s\n', savePath);
        if options.ExportSVG
            print(figFRET, [fullfile(outPath, [baseName '_FRET']) '.svg'], '-dsvg');
        end
    end

    %% --- Classification image -----------------------------------------
    if options.Class
        M = T.MaskClass{rowIdx};
        cats = unique(M(M > 0));   % categories present, excluding background
        nCats = numel(cats);

        % Build a colormap: row 1 = background (black), rows 2..end = category colours
        maxCat  = max(cats);
        cmap    = zeros(maxCat + 1, 3);          % row 1 = bg (black)
        for k = 1:maxCat
            cmap(k + 1, :) = getColorClass(k);  % row k+1 = category k
        end

        figClass = figure('Name', [baseName ' – Classification']);
        imagesc(M, [0 maxCat]);
        colormap(figClass, cmap);
        axis square;
        axis off;
        title(strrep(baseName, '_', '\_'), 'FontWeight', 'normal');

        % Manual legend patches (colorbar not meaningful for categorical data)
        hold on;
        hPatches = gobjects(nCats, 1);
        for k = 1:nCats
            hPatches(k) = patch(NaN, NaN, getColorClass(cats(k)), ...
                                'EdgeColor', 'none');
        end
        legend(hPatches, arrayfun(@(k) sprintf('Cat%d', k), cats, ...
               'UniformOutput', false), ...
               'Location', 'southoutside', 'Orientation', 'horizontal', ...
               'TextColor', 'white', 'Color', 'none', 'EdgeColor', 'none');
        hold off;

        % save
        savePath = fullfile(outPath, [baseName '_Class.tif']);
        exportgraphics(figClass, savePath, 'Resolution', 300);
        fprintf('Saved Class image: %s\n', savePath);
        if options.ExportSVG
            print(figClass, [fullfile(outPath, [baseName '_Class']) '.svg'], '-dsvg');
            % exportgraphics(figClass, [fullfile(outPath, [baseName '_Class']) '.pdf'], 'ContentType', 'vector');
        end
    end

        %% --- Skeleton image -----------------------------------------
    if options.Skel

        S = T.ThicknessMapSkel{rowIdx};

        % dilate to 3-pixel thickness
        se      = strel('disk', 1);          % radius 1  →  ~3 px wide
        skel3px = imdilate(S, se);

        % replace zeros with NaN
        skel_nan          = double(skel3px);
        skel_nan(skel_nan == 0) = NaN;

        % figure
        figSkel = figure('Name', [baseName ' – Skel']);
        imagesc(skel_nan, 'AlphaData', ~isnan(skel_nan)); % makes NaN pixels transparant
        % set colormap
        % Add black as the lowest color in the map
        cmap = cool(256);
        cmap = [0 0 0; cmap];          % prepend black
        % Set NaN pixels to a value just below the data minimum
        skel_display = skel_nan;
        skel_display(isnan(skel_nan)) = 0;   % 0 will map to black
        colormap(figSkel, cmap);
        imagesc(skel_display);
        fprintf('Colorbar max before ThickMax override: %.4f\n', max(skel_display(:)));
        clim([0 , options.ThickMax]);    % 0 → black, data range → cool
        % set(gca, 'Color', 'k');        % NaN pixels → black
        % set(gcf, 'Color', 'k');        % figure background → black too (optional)
        % colormap(figSkel, 'autumn');
        c = colorbar;
        c.Label.String = 'Thickness (µm)';
        % c.Color = 'w'; % if surrounding image bg is black
        axis square;
        axis off;
        title(strrep(baseName, '_', '\_'), 'FontWeight', 'normal');

        % save
        savePath = fullfile(outPath, [baseName '_Skel.tif']);
        exportgraphics(figSkel, savePath, 'Resolution', 300);
        fprintf('Saved Skel image: %s\n', savePath);
        if options.ExportSVG
            print(figSkel, [fullfile(outPath, [baseName '_Skel']) '.svg'], '-dsvg');
            % exportgraphics(figClass, [fullfile(outPath, [baseName '_Class']) '.pdf'], 'ContentType', 'vector');
        end
    end

end % main function

%% =========================================================================
function BlackJet = loadBlackJet()
% Loads BlackJet colormap. Expects BlackJet.mat on the MATLAB path.
    S = load('BlackJet.mat');
    fields = fieldnames(S);
    BlackJet = S.(fields{1});   % grab whatever variable is stored in the .mat
end


% function PaperImg_ExportIntFRETimage(T, Input, SampleFolder, EntryNumber, options)
% % ExportIntFRETimage  Export Int and/or FRET images as .tif for a given
% %                    sample/entry combination found in table T.
% %
% % Usage:
% %   ExportIntFRETimage(T, Input, SampleFolder, EntryNumber)
% %   ExportIntFRETimage(T, Input, SampleFolder, EntryNumber, options)
% %
% % Inputs:
% %   T             - Table containing at least columns: SampleFolder,
% %                   EntryNumber, IntD, FRETimageFA (cell arrays of NxN doubles)
% %   Input         - Struct with fields: OutPath, PixSize, Image_Resolution
% %   SampleFolder  - String/char, value to match in T.SampleFolder
% %   EntryNumber   - Scalar, value to match in T.EntryNumber
% %
% % Optional name-value inputs (options):
% %   options.Int          - Logical, export Int image  (default: true)
% %   options.FRET         - Logical, export FRET image (default: true)
% %   options.Class        - Logical, export Class image (default: false)
% %   options.Skel         - Logical, export Skel image (default: false)
% %   options.LegendInside - Logical, overlay colorbar/legend ON the image
% %                          matrix in data coordinates (default: true).
% %                          Scale bar black bg is always applied.
% 
%     arguments
%         T               table
%         Input           struct
%         SampleFolder    {mustBeTextScalar}
%         EntryNumber     (1,1) double
%         options.Int          (1,1) logical = true
%         options.FRET         (1,1) logical = true
%         options.Class        (1,1) logical = false
%         options.Skel         (1,1) logical = false
%         options.LegendInside (1,1) logical = true
%         options.IntLow       (1,1) double  = 1
%         options.IntHigh      (1,1) double  = 99
%         options.ThickMax     (1,1) double  = inf
%         options.ScaleBarUm   (1,1) double  = 20
%         options.ExportSVG    (1,1) logical = false
%     end
% 
%     %% --- Locate the unique row in T -----------------------------------
%     rowMask = strcmp(T.SampleFolder, SampleFolder) & ...
%               (T.EntryNumber == EntryNumber);
% 
%     nHits = sum(rowMask);
%     if nHits == 0
%         error('ExportIntFRETimage:noMatch', ...
%               'No row found for SampleFolder="%s" and EntryNumber=%d.', ...
%               SampleFolder, EntryNumber);
%     elseif nHits > 1
%         error('ExportIntFRETimage:ambiguousMatch', ...
%               'Multiple rows (%d) found for SampleFolder="%s" and EntryNumber=%d.', ...
%               nHits, SampleFolder, EntryNumber);
%     end
% 
%     rowIdx   = find(rowMask);
%     outPath  = Input.OutPath;
%     baseName = sprintf('%s_Cell_%d', SampleFolder, EntryNumber);
%     dim      = double(Input.Image_Resolution);   % e.g. 512
% 
%     %% --- Intensity image ----------------------------------------------
%     if options.Int
%         A    = T.IntD{rowIdx};
%         low  = double(prctile(A(:), options.IntLow));
%         high = double(prctile(A(:), options.IntHigh));
% 
%         figInt = figure('Name', [baseName ' – Intensity']);
%         ax = axes('Parent', figInt);
%         imagesc(ax, A, [low high]);
%         colormap(figInt, gray);
%         axis(ax, 'square');
%         axis(ax, 'off');
%         title(ax, strrep(baseName, '_', '\_'), 'FontWeight', 'normal');
% 
%         % Scale bar (always with black background box)
%         addScaleBar(ax, Input, options.ScaleBarUm, dim);
% 
%         % Colorbar overlaid on the image matrix
%         if options.LegendInside
%             addInsideColorbar(ax, gray(256), [low high], '', dim);
%         else
%             colorbar(ax);
%         end
% 
%         savePath = fullfile(outPath, [baseName '_Int.tif']);
%         exportgraphics(figInt, savePath, 'Resolution', 300);
%         fprintf('Saved Int image: %s\n', savePath);
%         if options.ExportSVG
%             print(figInt, fullfile(outPath, [baseName '_Int.svg']), '-dsvg');
%         end
%     end
% 
%     %% --- FRET image ---------------------------------------------------
%     if options.FRET
%         BlackJet = loadBlackJet();
%         F = T.FRETimageFA{rowIdx};
% 
%         figFRET = figure('Name', [baseName ' – FRET']);
%         ax = axes('Parent', figFRET);
%         imagesc(ax, F, [0 40]);
%         colormap(figFRET, BlackJet);
%         axis(ax, 'square');
%         axis(ax, 'off');
%         title(ax, strrep(baseName, '_', '\_'), 'FontWeight', 'normal');
% 
%         if options.LegendInside
%             addInsideColorbar(ax, BlackJet, [0 40], 'Apparent FRET %', dim);
%         else
%             c = colorbar(ax);
%             c.Label.String = 'Apparent FRET Efficiency %';
%         end
% 
%         savePath = fullfile(outPath, [baseName '_FRET.tif']);
%         exportgraphics(figFRET, savePath, 'Resolution', 300);
%         fprintf('Saved FRET image: %s\n', savePath);
%         if options.ExportSVG
%             print(figFRET, fullfile(outPath, [baseName '_FRET.svg']), '-dsvg');
%         end
%     end
% 
%     %% --- Classification image -----------------------------------------
%     if options.Class
%         M      = T.MaskClass{rowIdx};
%         cats   = unique(M(M > 0));
%         nCats  = numel(cats);
%         maxCat = max(cats);
% 
%         cmap = zeros(maxCat + 1, 3);
%         for k = 1:maxCat
%             cmap(k + 1, :) = getColorClass(k);
%         end
% 
%         figClass = figure('Name', [baseName ' – Classification']);
%         ax = axes('Parent', figClass);
%         imagesc(ax, M, [0 maxCat]);
%         colormap(figClass, cmap);
%         axis(ax, 'square');
%         axis(ax, 'off');
%         title(ax, strrep(baseName, '_', '\_'), 'FontWeight', 'normal');
% 
%         if options.LegendInside
%             catNames  = arrayfun(@(k) sprintf('Cat%d', k), cats, ...
%                                  'UniformOutput', false);
%             catColors = cell2mat(arrayfun(@(k) getColorClass(k), cats, ...
%                                          'UniformOutput', false)');
%             addInsideLegend(ax, catNames, catColors, dim);
%         else
%             hold(ax, 'on');
%             hPatches = gobjects(nCats, 1);
%             for k = 1:nCats
%                 hPatches(k) = patch(ax, NaN, NaN, getColorClass(cats(k)), ...
%                                     'EdgeColor', 'none');
%             end
%             legend(ax, hPatches, ...
%                    arrayfun(@(k) sprintf('Cat%d', k), cats, ...
%                             'UniformOutput', false), ...
%                    'Location', 'southoutside', 'Orientation', 'horizontal', ...
%                    'TextColor', 'white', 'Color', 'none', 'EdgeColor', 'none');
%             hold(ax, 'off');
%         end
% 
%         savePath = fullfile(outPath, [baseName '_Class.tif']);
%         exportgraphics(figClass, savePath, 'Resolution', 300);
%         fprintf('Saved Class image: %s\n', savePath);
%         if options.ExportSVG
%             print(figClass, fullfile(outPath, [baseName '_Class.svg']), '-dsvg');
%         end
%     end
% 
%     %% --- Skeleton image -----------------------------------------------
%     if options.Skel
%         S = T.ThicknessMapSkel{rowIdx};
% 
%         se           = strel('disk', 1);
%         skel3px      = imdilate(S, se);
%         skel_plot    = double(skel3px);
%         skel_plot(skel_plot == 0) = NaN;
% 
%         cmap = cool(256);
%         cmap = [0 0 0; cmap];
% 
%         skel_draw           = skel_plot;
%         skel_draw(isnan(skel_plot)) = 0;
% 
%         figSkel = figure('Name', [baseName ' – Skel']);
%         ax = axes('Parent', figSkel);
%         imagesc(ax, skel_draw);
%         colormap(figSkel, cmap);
%         clim(ax, [0, options.ThickMax]);
%         axis(ax, 'square');
%         axis(ax, 'off');
%         title(ax, strrep(baseName, '_', '\_'), 'FontWeight', 'normal');
% 
%         if options.LegendInside
%             addInsideColorbar(ax, cmap, [0, double(options.ThickMax)], ...
%                               'Thickness (µm)', dim);
%         else
%             c = colorbar(ax);
%             c.Label.String = 'Thickness (µm)';
%         end
% 
%         savePath = fullfile(outPath, [baseName '_Skel.tif']);
%         exportgraphics(figSkel, savePath, 'Resolution', 300);
%         fprintf('Saved Skel image: %s\n', savePath);
%         if options.ExportSVG
%             print(figSkel, fullfile(outPath, [baseName '_Skel.svg']), '-dsvg');
%         end
%     end
% 
% end % main function
% 
% 
% %%%========================================================================
% %%%  LOCAL HELPER FUNCTIONS  (all work in image data coordinates 1..dim)
% %%%========================================================================
% 
% %% -----------------------------------------------------------------------
% function addScaleBar(ax, Input, scaleBarUm, dim)
% % Draws a scale bar with a black background box directly on the image axes.
% % Everything is in data (pixel) coordinates: x in [1,dim], y in [1,dim].
% % Draw order: black bg first -> white bar on top -> white text on top.
% 
%     pixSize     = double(Input.PixSize);
%     barLengthPx = scaleBarUm / pixSize;
% 
%     marginRight  = 22;   % px from right edge of image
%     marginBottom = 22;   % px from bottom edge of image
%     barThick     = 5;    % bar height in pixels
%     pad          = 7;    % padding around combined bar+text block
% 
%     % Bar rectangle (bottom-right, fully inset)
%     barX = dim - marginRight - barLengthPx;   % left edge of bar
%     barY = dim - marginBottom - barThick;     % top edge of bar
% 
%     % Text: centred above bar, a few pixels of gap
%     txtGap = 3;    % gap between top of bar and bottom of text
%     txtX   = barX + barLengthPx / 2;
%     txtY   = barY - txtGap;    % VerticalAlignment 'bottom' -> text sits above this
% 
%     % --- Estimate text height in data units ---
%     % Draw text first (off-screen) to get Extent, then reposition bg
%     hold(ax, 'on');
% 
%     hTxt = text(ax, txtX, txtY, sprintf('%g µm', scaleBarUm), ...
%         'Color', 'white', 'FontSize', 9, 'FontWeight', 'bold', ...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
%         'Units', 'data');
% 
%     drawnow;   % force renderer so Extent is populated
%     ext  = get(hTxt, 'Extent');   % [left bottom width height] in data coords
%     % ext(2) = bottom of text bounding box, ext(4) = height
% 
%     % Bounding box enclosing both text and bar
%     allLeft   = min(ext(1),            barX);
%     allBottom = min(ext(2),            barY);
%     allRight  = max(ext(1) + ext(3),   barX + barLengthPx);
%     allTop    = max(ext(2) + ext(4),   barY + barThick);
% 
%     % Add padding, clamp to image boundaries
%     bgLeft   = max(allLeft   - pad, 1);
%     bgBottom = max(allBottom - pad, 1);
%     bgRight  = min(allRight  + pad, dim);
%     bgTop    = min(allTop    + pad, dim);
% 
%     % 1) Black background (drawn first = underneath)
%     fill(ax, ...
%         [bgLeft bgRight bgRight bgLeft], ...
%         [bgBottom bgBottom bgTop bgTop], ...
%         'black', 'EdgeColor', 'none');
% 
%     % 2) White bar (on top of bg)
%     fill(ax, ...
%         [barX, barX+barLengthPx, barX+barLengthPx, barX], ...
%         [barY, barY, barY+barThick, barY+barThick], ...
%         'white', 'EdgeColor', 'none');
% 
%     % 3) Re-draw text so it sits above both bg and bar
%     % (delete the measurement copy and draw fresh)
%     delete(hTxt);
%     text(ax, txtX, txtY, sprintf('%g µm', scaleBarUm), ...
%         'Color', 'white', 'FontSize', 9, 'FontWeight', 'bold', ...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
%         'Units', 'data');
% 
%     hold(ax, 'off');
% end
% 
% 
% %% -----------------------------------------------------------------------
% function addInsideColorbar(ax, cmap, clims, labelStr, dim)
% % Overlays a vertical colorbar on the image axes using data coordinates.
% % Position: top-right corner of the image matrix.
% % Height = 25% of image, width ~4% of image.
% %
% % Strategy: draw a black filled rectangle first, then render the colormap
% % as a true-colour image patch on top, then draw tick lines + labels.
% 
%     nColors  = size(cmap, 1);
% 
%     % --- Geometry in data (pixel) coordinates -------------------------
%     marginRight = 18;
%     marginTop   = 18;
%     cbW         = round(dim * 0.04);     % ~4 % of image width
%     cbH         = round(dim * 0.25);     % 25 % of image height
%     pad         = 6;                     % padding around cb inside black box
%     tickLen     = 4;                     % tick mark length in pixels
%     fontSize    = 7;
% 
%     % Colorbar top-right position (image coords: y increases downward)
%     cbRight  = dim - marginRight;
%     cbLeft   = cbRight - cbW;
%     cbTop    = marginTop;
%     cbBottom = cbTop + cbH;
% 
%     % Black background box (includes space for ticks and label to the left)
%     labelW  = round(dim * 0.07);   % estimate for tick labels width
%     bgLeft   = cbLeft  - tickLen - labelW - pad;
%     bgRight  = cbRight + pad;
%     bgTop2   = cbTop   - pad - (fontSize + 4);  % headroom for title
%     bgBottom2= cbBottom + pad;
% 
%     % Clamp to image
%     bgLeft   = max(bgLeft,   1);
%     bgTop2   = max(bgTop2,   1);
%     bgRight  = min(bgRight,  dim);
%     bgBottom2= min(bgBottom2,dim);
% 
%     hold(ax, 'on');
% 
%     % 1) Black background
%     fill(ax, ...
%         [bgLeft bgRight bgRight bgLeft], ...
%         [bgTop2 bgTop2 bgBottom2 bgBottom2], ...
%         'black', 'EdgeColor', 'none');
% 
%     % 2) Colormap strip as individual horizontal filled rectangles
%     %    (avoids needing a nested axes or image call)
%     rowH = cbH / nColors;
%     for k = 1:nColors
%         y1 = cbTop    + (k-1) * rowH;
%         y2 = cbTop    + k     * rowH;
%         fill(ax, ...
%             [cbLeft cbRight cbRight cbLeft], ...
%             [y1     y1      y2      y2], ...
%             cmap(k,:), 'EdgeColor', 'none');
%     end
% 
%     % 3) Border around colorbar strip
%     plot(ax, [cbLeft cbRight cbRight cbLeft cbLeft], ...
%              [cbTop  cbTop   cbBottom cbBottom cbTop], ...
%         'w-', 'LineWidth', 0.5);
% 
%     % 4) Ticks and labels (on the LEFT side of the strip)
%     nTicks   = 5;
%     tickVals = linspace(double(clims(1)), double(clims(2)), nTicks);
%     tickYpos = linspace(cbTop, cbBottom, nTicks);
% 
%     for t = 1:nTicks
%         % Tick line
%         plot(ax, [cbLeft - tickLen, cbLeft], ...
%                  [tickYpos(t), tickYpos(t)], ...
%             'w-', 'LineWidth', 0.5);
%         % Tick label
%         text(ax, cbLeft - tickLen - 2, tickYpos(t), ...
%             sprintf('%.4g', tickVals(t)), ...
%             'Color', 'white', 'FontSize', fontSize, ...
%             'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
%     end
% 
%     % 5) Optional title above the colorbar
%     if ~isempty(labelStr)
%         text(ax, (cbLeft + cbRight)/2, bgTop2 + 2, labelStr, ...
%             'Color', 'white', 'FontSize', fontSize, ...
%             'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
%             'Interpreter', 'none');
%     end
% 
%     hold(ax, 'off');
% end
% 
% 
% %% -----------------------------------------------------------------------
% function addInsideLegend(ax, catNames, catColors, dim)
% % Overlays a categorical legend on the image axes in data coordinates.
% % Position: top-right corner, vertically organised, black background.
% 
%     nCats    = numel(catNames);
%     fontSize = 8;
% 
%     % --- Geometry ---
%     marginRight = 18;
%     marginTop   = 18;
%     swatchW     = round(dim * 0.04);
%     swatchH     = round(dim * 0.03);
%     rowGap      = round(dim * 0.01);
%     textGap     = 5;
%     pad         = 6;
% 
%     totalH = nCats * swatchH + (nCats-1) * rowGap;
% 
%     bgRight  = dim - marginRight;
%     bgTop    = marginTop;
% 
%     hold(ax, 'on');
% 
%     % Estimate text width for background sizing (rough: ~6 px per char at fontSize 8)
%     maxChars = max(cellfun(@numel, catNames));
%     textW    = maxChars * 6 + textGap;
% 
%     bgLeft   = bgRight - swatchW - textW - 2*pad;
%     bgBottom = bgTop   + totalH  + 2*pad;
% 
%     % Clamp
%     bgLeft   = max(bgLeft,   1);
%     bgBottom = min(bgBottom, dim);
% 
%     % 1) Black background
%     patch(ax, ...
%         [bgLeft,  bgRight, bgRight, bgLeft], ...
%         [bgTop,   bgTop,   bgBottom, bgBottom], ...
%         'k', 'EdgeColor', 'none');
% 
%     % 2) Swatches + labels
%     for k = 1:nCats
%         rowTop = bgTop + pad + (k-1) * (swatchH + rowGap);
% 
%         swLeft = bgLeft + pad;
%         swRight= swLeft + swatchW;
% 
%         % Colour swatch
%         patch(ax, ...
%             [swLeft,  swRight, swRight, swLeft], ...
%             [rowTop,  rowTop,  rowTop+swatchH,  rowTop+swatchH], ...
%             'w', ...                          % placeholder color (overridden below)
%             'FaceColor', catColors(k,:), ...
%             'EdgeColor', 'none');
% 
%         % Label
%         text(ax, swRight + textGap, rowTop + swatchH/2, catNames{k}, ...
%             'Color', 'white', 'FontSize', fontSize, ...
%             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
%             'Interpreter', 'none');
%     end
% 
%     hold(ax, 'off');
% end
% 
% 
% %% -----------------------------------------------------------------------
% function BlackJet = loadBlackJet()
%     S = load('BlackJet.mat');
%     fields = fieldnames(S);
%     BlackJet = S.(fields{1});
% end
