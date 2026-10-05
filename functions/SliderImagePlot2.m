function [finalValue] = SliderImagePlot2(Int)
% SliderImagePlot lets you adjust threshold by slider or typing a value

    Int = double(Int);
    a = 1; % initial value

    fig = figure('Position', [100, 100, 800, 600]);

    % Set the initial value of 'a' as the application data
    setappdata(fig, 'FinalValue', a);

    % Create an axes for the imagesc plot
    ax = axes('Parent', fig, 'Position', [0.1, 0.3, 0.8, 0.6]);
    
    % Create the initial imagesc plot
    hImage = imagesc(ax, imbinarize((Int/max(Int(:))),(a/max(Int(:)))));
  axis(ax, 'image', 'off');
    colorbar;

    % Slider setup
    maxi = max(Int(:));
    slider = uicontrol('Style', 'slider', 'Parent', fig, ...
        'Position', [75, 50, 500, 20], 'Min', 1, 'Max', maxi, ...
        'Value', a, 'SliderStep', [0.01, 0.1]);

    % Numeric edit box next to slider
    editBox = uicontrol('Style','edit','Parent',fig,...
        'Position',[600, 50, 80, 25],...
        'String',num2str(a),...
        'BackgroundColor','white',...
        'Callback',@editBoxCallback);

    % Button to close
    closeButton = uicontrol('Style', 'pushbutton', 'Parent', fig, ...
        'Position', [350, 10, 100, 30], 'String', 'Close & Save', ...
        'Callback', @closeAndSave);

    % Slider labels
    LabelMin = uicontrol('Parent',fig,'Style','text','Position',[30,54,23,23],...
                'String','0');
    LabelMax = uicontrol('Parent',fig,'Style','text','Position',[560,54,40,23],...
                'String',num2str(maxi));

    % Add a callback function to the slider
    addlistener(slider, 'Value', 'PostSet', @(src, event) updatePlot(round(get(slider,'Value'))));

    % Nested update function
    function updatePlot(val)
        a = round(val);
        set(hImage, 'CData', imbinarize((Int/max(Int(:))),(a/max(Int(:)))));
        title(ax, sprintf('Intensity Threshold = %d', a));
        % Sync edit box with slider
        set(editBox,'String',num2str(a));
        drawnow;
    end

    % Callback when user types in the edit box
    function editBoxCallback(src,~)
        val = str2double(get(src,'String'));
        if isnan(val) || val < 1 || val > maxi
            % Reset if invalid
            set(src,'String',num2str(round(get(slider,'Value'))));
        else
            set(slider,'Value',val);
            updatePlot(val);
        end
    end

    % Callback for the button
    function closeAndSave(~, ~)
        finalValue = round(get(slider, 'Value'));
        setappdata(fig, 'FinalValue', finalValue);
        uiresume(fig);   % just resume, don't close
    end

    % --- Wait for user and then return value ---
    uiwait(fig);  
    finalValue = getappdata(fig,'FinalValue');
    if ishandle(fig)
        close(fig);
    end

end
