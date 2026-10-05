function Plot_ColorFAMask(T, SampleFolder, Nr)
% Plot_ColorFAMask
% --------------------------------------------------------
% Visualizes T.MaskClass for a specific table entry.
%
% INPUTS:
%   T            = main table
%   SampleFolder = string/pattern to match T.SampleFolder
%   Nr           = EntryNumber (unique per sample)
%
% Colors:
%   0 = black (background)
%   1 = blue    (Cat1)
%   2 = green   (Cat2)
%   3 = orange  (Cat3)
%   4 = purple  (Cat4)
% --------------------------------------------------------

    %% --- Find the correct row in T ---
    mask = contains(T.SampleFolder, SampleFolder) & T.EntryNumber == Nr;
    
    if ~any(mask)
        error('No matching row found for SampleFolder "%s" and EntryNumber %d.', ...
            SampleFolder, Nr);
    end
    
    rowID = find(mask, 1, 'first');


    
    %% --- Get classification mask ---
    ClassMask = T.MaskClass{rowID};   % 512×512 double with values 0–4
    
    if isempty(ClassMask)
        error('T.MaskClass{%d} is empty.', rowID);
    end
    
    %% --- Build a colormap (0–4) ---
    % IMPORTANT: colormap rows correspond to integer category values
    % index 1 → value 0
    % index 2 → value 1
    % etc.
    % Build colormap using getColorClass()
    cmap = [
        0 0 0;                           % 0 = background
        getColorClass('Cat1');           % 1 = Cat1
        getColorClass('Cat2');           % 2 = Cat2
        getColorClass('Cat3');           % 3 = Cat3
        getColorClass('Cat4')            % 4 = Cat4
    ];
    
    %% --- Plot THREE images side by side ---
    figure;
    tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');
    
    
    % -------- Panel 1: T.IntD ---------
    nexttile;
    if ismember('IntD', T.Properties.VariableNames)
        I = T.IntD{rowID};
    else
        I = T.Int{rowID};
    end
    hINT = imagesc(I);
    axis image off;
    title('Intensity (IntD)', 'Interpreter','none');
    colormap(gca, 'gray');
    colorbar;
    
    % Apply your CLim-based brightness control
    ax = gca;
    clim = ax.CLim;
    origMin = clim(1);     % preserve original black point
    autoMax = 100;         % match your other function
    set(ax, 'CLim', [origMin autoMax]);
    
    % -------- Panel 2: T.FRETimageFA ---------
    load('BlackJet.mat');
    nexttile;
    FRET = T.FRETimageFA{rowID};
    Fdisp = FRET;
    Fdisp(~isfinite(Fdisp)) = 0;   % NaN background -> 0
    
    imagesc(Fdisp);
    colormap(gca, BlackJet);        % your custom colormap
    colorbar;
    
    % upper limit based on nonzero values
    maxF = 40;
    if isempty(maxF) || maxF <= 0 || ~isfinite(maxF)
        maxF = 1;   % fallback
    end
    
    caxis([0 maxF]);                % same style as your FRET visualization
    title('FRET Image (FRETimageFA)');
    axis image off;
    
    
    % -------- Panel 3: T.MaskClass ---------
    nexttile;
    imagesc(ClassMask);
    axis image off;
    colormap(gca, cmap);         % <-- using getColorClass-based cmap
    caxis([0 4]);
    title('FA Classification (MaskClass)');
    % --- Colored legend using getColorClass() ---
    hold on;
    hCat1 = plot(nan, nan, 's', 'MarkerFaceColor', getColorClass('Cat1'),  'MarkerEdgeColor','none');
    hCat2 = plot(nan, nan, 's', 'MarkerFaceColor', getColorClass('Cat2'),  'MarkerEdgeColor','none');
    hCat3 = plot(nan, nan, 's', 'MarkerFaceColor', getColorClass('Cat3'),  'MarkerEdgeColor','none');
    hCat4 = plot(nan, nan, 's', 'MarkerFaceColor', getColorClass('Cat4'),  'MarkerEdgeColor','none');
    hold off;
    legend([hCat1 hCat2 hCat3 hCat4], ...
        {'Cat1','Cat2','Cat3','Cat4'}, ...
        'Location','southoutside', 'NumColumns',4);

end