function mask = interactive_roi(img, maskOverlay) % bg should be 0 (?) if not, please do this beforehand
% INTERACTIVE_ROI lets you manually select multiple ROIs and combine them.
% Displays both the original image and overlayed mask side by side.
%
% mask = INTERACTIVE_ROI(img, maskOverlay)
% img  - 2D numeric array (e.g. uint16 microscopy image)
% maskOverlay - optional logical array same size as img
%
% Returns:
% mask - logical array of selected ROIs

    if ~ismatrix(img)
        error('Input must be a 2D image.');
    end
    if nargin < 2
        maskOverlay = [];
    end

    mask = false(size(img));

    % Create figure
    fig = figure('Name','Draw ROIs — Click "Done" when finished', ...
                 'NumberTitle','off', 'Units','normalized', 'Position',[0.1 0.1 0.9 0.9]);

    % Create two axes: left = original, right = overlay
    gap = 0.02;      % gap between images (2% of figure width)
    leftWidth = (1 - 3*gap)/2;  % width of each image
    
    % Left image
    ax1 = axes('Parent', fig, 'Position', [gap, 0.1, leftWidth, 0.8]);
    imshow(img, [], 'Parent', ax1, 'InitialMagnification', 'fit');
    title(ax1, 'Original Image');
    
    % Right image
    ax2 = axes('Parent', fig, 'Position', [2*gap + leftWidth, 0.1, leftWidth, 0.8]);
    hImg = imshow(img, [], 'Parent', ax2, 'InitialMagnification', 'fit');
    title(ax2, 'Overlayed Image');
    hold(ax2, 'on');

    % Overlay mask if provided
    if ~isempty(maskOverlay)
        if ~isequal(size(maskOverlay), size(img))
            error('maskOverlay must be the same size as img.');
        end
        if ~islogical(maskOverlay)
            maskOverlay = mat2gray(maskOverlay);
        end
        % maskOverlay = ~maskOverlay;  % invert if necessary
        overlayRGB = cat(3, zeros(size(maskOverlay)), maskOverlay, zeros(size(maskOverlay)));
        hOverlay = imshow(overlayRGB, 'Parent', ax2);
        set(hOverlay, 'AlphaData', 0.4 * maskOverlay);
    end

    % --- Display scaling for both axes ---
    clim1 = get(ax1, 'CLim');
    origMin1 = clim1(1);
    clim2 = get(ax2, 'CLim');    
    origMin2 = clim2(1);
    maxSliderRange = double(max(img(:)));
    %autoMax = maxSliderRange/2; % autoMax = 200;
    autoMax = maxSliderRange/4; % for SarahV data
    set(ax1, 'CLim', [origMin1 autoMax]);
    set(ax2, 'CLim', [origMin2 autoMax]);
    % store info for slider
    setappdata(fig, 'maxSliderRange', maxSliderRange);
    setappdata(fig, 'stopDrawing', false);
    setappdata(fig, 'axLeft', ax1);
    setappdata(fig, 'axRight', ax2);
    setappdata(fig, 'origMinLeft', origMin1);
    setappdata(fig, 'origMinRight', origMin2);
    setappdata(fig, 'autoMax', autoMax);

    % %old: Left image (original)
    % clim1 = get(ax1, 'CLim');    
    % origMin1 = clim1(1);         
    % newMaxLeft = 50;            
    % set(ax1, 'CLim', [origMin1 newMaxLeft]);
    % % Right image (overlay)
    % clim2 = get(ax2, 'CLim');    
    % origMin2 = clim2(1);         
    % newMaxRight = 50;           
    % set(ax2, 'CLim', [origMin2 newMaxRight]);
    % 
    % % --- Store drawing state for ROI ---
    % setappdata(fig, 'stopDrawing', false);
    % % setappdata(fig, 'axHandleRight', ax1);
    % % setappdata(fig, 'origMinRight', origMin1);
    % % setappdata(fig, 'axHandleLeft', ax1);
    % % setappdata(fig, 'origMinLeft', origMin1);
    % % %setappdata(fig, 'origMin', origMin2);
    % setappdata(fig, 'axLeft', ax1);
    % setappdata(fig, 'axRight', ax2);
    % setappdata(fig, 'origMinLeft', origMin1);
    % setappdata(fig, 'origMinRight', origMin2);

    % --- "Done" button ---
    uicontrol('Style','pushbutton', ...
              'String','Done', ...
              'FontSize',10, ...
              'Units','normalized', ...
              'Position',[0.85 0.02 0.1 0.05], ...
              'Callback',@(src,evt) setappdata(fig, 'stopDrawing', true));

    % --- Contrast Slider (0–1) ---
    autoMax = getappdata(fig, 'autoMax');
    % Slider range: 0 → dark, 1 → full brightness
    % Starting value matches the autoMax we just applied
    %maxSliderRange = 750; % adjust if images have different intensity scales
    startValue = autoMax / maxSliderRange;
    uicontrol('Style', 'slider', ...
              'Min', 0, 'Max', 1, 'Value', startValue, ...
              'Units', 'normalized', ...
              'Position', [0.35 0.02 0.22 0.05], ...
              'Callback', @brightnessSliderCallback);
    
    uicontrol('Style', 'text', ...
              'String', 'Brightness', ...
              'Units', 'normalized', ...
              'Position', [0.35 0.07 0.22 0.03], ...
              'HorizontalAlignment', 'center');
    % % --- "Set Max" buttons ---
    % uicontrol('Style','pushbutton', ...
    %           'String','Set Max (Right)', ...
    %           'FontSize',10, ...
    %           'Units','normalized', ...
    %           'Position',[0.72 0.02 0.1 0.05], ...
    %           'Callback',@setMaxRightCallback);
    % 
    % uicontrol('Style','pushbutton', ...
    %       'String','Set Max (Left)', ...
    %       'FontSize',10, ...
    %       'Units','normalized', ...
    %       'Position',[0.59 0.02 0.12 0.05], ...
    %       'Callback',@setMaxLeftCallback);

    % Main drawing loop
    while ishandle(fig) && ~getappdata(fig, 'stopDrawing')
        h = drawfreehand(ax2, 'Color', 'r');  % draw on overlay
        if ~ishandle(h) || getappdata(fig, 'stopDrawing')
            break;
        end

        roiMask = createMask(h, hImg);
        mask = mask | roiMask;

        % show ROI boundary
        B = bwboundaries(roiMask, 'noholes');
        for k = 1:numel(B)
            b = B{k};
            plot(ax2, b(:,2), b(:,1), 'g', 'LineWidth', 1);
        end
        drawnow;
    end

    if ishandle(fig)
        close(fig);
    end
end

% --- Button callback ---
function brightnessSliderCallback(src, ~)
    fig = ancestor(src, 'figure');

    axL = getappdata(fig, 'axLeft');
    axR = getappdata(fig, 'axRight');

    origMinL = getappdata(fig, 'origMinLeft');
    origMinR = getappdata(fig, 'origMinRight');

    % same scale as the auto brightness
    % maxSliderRange = 750;
    maxSliderRange = getappdata(fig, 'maxSliderRange');

    sliderValue = get(src, 'Value');
    newMax = sliderValue * maxSliderRange;
    
    if newMax < origMinL
        newMax = origMinL + eps;
    end

    set(axL, 'CLim', [origMinL newMax]);
    set(axR, 'CLim', [origMinR newMax]);
end
% function setMaxRightCallback(~, ~)
%     fig = gcbf; % current figure
%     ax = getappdata(fig, 'axHandleRight');
%     origMin = getappdata(fig, 'origMinRight');
% 
%     answer = inputdlg('Enter new display maximum value:', 'Set Display Max (Right)', [1 40]);
%     if isempty(answer)
%         return;
%     end
% 
%     newMax = str2double(answer{1});
%     if isnan(newMax) || newMax <= origMin
%         warndlg('Invalid maximum value.');
%         return;
%     end
% 
%     set(ax, 'CLim', [origMin newMax]);
% end
% 
% function setMaxLeftCallback(~, ~)
%     fig = gcbf;
%     ax = getappdata(fig, 'axHandleLeft');
%     origMin = getappdata(fig, 'origMinLeft');
% 
%     answer = inputdlg('Enter new display maximum value for LEFT image:', ...
%                       'Set Display Max (Left)', [1 40]);
%     if isempty(answer)
%         return;
%     end
% 
%     newMax = str2double(answer{1});
%     if isnan(newMax) || newMax <= origMin
%         warndlg('Invalid maximum value.');
%         return;
%     end
% 
%     set(ax, 'CLim', [origMin newMax]);
% end

% function mask = interactive_roi(img, maskOverlay)
% % DRAW_COMBINED_MASK lets you manually select multiple ROIs and combine them.
% %   mask = DRAW_COMBINED_MASK(img)
% %   img  - 2D numeric array (e.g. uint16 microscopy image)
% %   mask - logical array same size as img, with selected regions = true
% %
% % Buttons:
% %   - "Done": finish ROI drawing and return mask (doesn't work right now)
% %   - "Set Max": change display scaling maximum for better visibility
% 
%     if ~ismatrix(img)
%         error('Input must be a 2D image.');
%     end
% 
%     if nargin < 2
%         maskOverlay = [];
%     end
% 
%     mask = false(size(img));
% 
%     % Create figure & axes
%     fig = figure('Name','Draw ROIs — Click "Done" when finished', ...
%                  'NumberTitle','off');
%     ax = axes('Parent',fig);
%     hImg = imshow(img, [], 'Parent', ax, 'InitialMagnification', 'fit');
%     title(ax, 'Draw ROIs.');
%     hold(ax, 'on');
% 
%     % If maskOverlay is provided, display it semi-transparent
%     if ~isempty(maskOverlay)
%         % Ensure maskOverlay matches img size
%         if ~isequal(size(maskOverlay), size(img))
%             error('maskOverlay must be the same size as img.');
%         end
% 
%         % Normalize to [0,1] if not logical
%         if ~islogical(maskOverlay)
%             maskOverlay = mat2gray(maskOverlay);
%         end
% 
%         % Invert maskOverlay
%         maskOverlay = ~maskOverlay;
% 
%         % Create a green overlay with transparency
%         overlayRGB = cat(3, zeros(size(maskOverlay)), maskOverlay, zeros(size(maskOverlay)));
%         hOverlay = imshow(overlayRGB, 'Parent', ax);
%         set(hOverlay, 'AlphaData', 0.4 * maskOverlay);  % adjust transparency as needed
%         % Create a red overlay with transparency
%         %overlayRGB = cat(3, maskOverlay, zeros(size(maskOverlay)), zeros(size(maskOverlay)));
%         %hOverlay = imshow(overlayRGB, 'Parent', ax);
%         %set(hOverlay, 'AlphaData', 0.8 * maskOverlay); % adjust transparency as desired
%     end
% 
%     % Store original scaling
%     clim = get(ax, 'CLim'); % [min max]
%     origMin = clim(1);
%     set(ax, 'CLim', [origMin 200]);   % set display max to 100 already
%     setappdata(fig, 'stopDrawing', false);
%     setappdata(fig, 'axHandle', ax);
%     setappdata(fig, 'origMin', origMin);
% 
%     % --- "Done" button ---
%     uicontrol('Style','pushbutton', ...
%               'String','Done', ...
%               'FontSize',10, ...
%               'Units','normalized', ...
%               'Position',[0.85 0.02 0.1 0.05], ...
%               'Callback',@(src,evt) setappdata(fig, 'stopDrawing', true));
% 
%     % --- "Set Max" button ---
%     uicontrol('Style','pushbutton', ...
%               'String','Set Max', ...
%               'FontSize',10, ...
%               'Units','normalized', ...
%               'Position',[0.72 0.02 0.1 0.05], ...
%               'Callback',@setMaxCallback);
% 
%     % Main drawing loop
%     while ishandle(fig) && ~getappdata(fig, 'stopDrawing')
%         h = drawfreehand(ax, 'Color', 'r');  % waits for user to finish ROI
%         if ~ishandle(h) || getappdata(fig, 'stopDrawing')
%             break;
%         end
% 
%         % Create ROI mask on the base image only
%         roiMask = createMask(h, hImg);
%         mask = mask | roiMask;
% 
%         % show ROI boundary
%         B = bwboundaries(roiMask, 'noholes');
%         for k = 1:numel(B)
%             b = B{k};
%             plot(ax, b(:,2), b(:,1), 'g', 'LineWidth', 1);
%         end
%         drawnow;
%     end
% 
%     if ishandle(fig)
%         close(fig);
%     end
% end
% 
% % --- Button callback to change brightness ---
% function setMaxCallback(~, ~)
%     fig = gcbf; % current figure
%     ax = getappdata(fig, 'axHandle');
%     origMin = getappdata(fig, 'origMin');
% 
%     answer = inputdlg('Enter new display maximum value:', ...
%                       'Set Display Max', [1 40]);
%     if isempty(answer)
%         return;
%     end
% 
%     newMax = str2double(answer{1});
%     if isnan(newMax) || newMax <= origMin
%         warndlg('Invalid maximum value.');
%         return;
%     end
% 
%     set(ax, 'CLim', [origMin newMax]);
% end