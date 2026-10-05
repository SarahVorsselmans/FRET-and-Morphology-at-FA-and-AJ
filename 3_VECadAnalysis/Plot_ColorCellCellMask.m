function Plot_ColorCellCellMask(T, SampleFolder, Nr)
% INPUTS:
%   T            = main table
%   SampleFolder = string/pattern to match T.SampleFolder
%   Nr           = EntryNumber (unique per sample)
% --------------------------------------------------------

    %% --- Find the correct row in T ---
    mask = contains(T.SampleFolder, SampleFolder) & T.EntryNumber == Nr;
    
    if ~any(mask)
        error('No matching row found for SampleFolder "%s" and EntryNumber %d.', ...
            SampleFolder, Nr);
    end
    
    rowID = find(mask, 1, 'first');
    
    %% --- Plot THREE images side by side ---
    figure;
    tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');
    
    
    % -------- Panel 1: T.IntD ---------
    nexttile;
    I = T.IntD{rowID};
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
    hasFRETimage = ismember('FRETimageFA', T.Properties.VariableNames) && ~isempty(T.FRETimageFA{rowID});
    hasFRETindex = ismember('FRETindex',   T.Properties.VariableNames) && ~isempty(T.FRETindex{rowID});

    if hasFRETimage
        % User-supplied binary mask (0/1)
        %FRET = ~isnan(T.FRETimageFA{rowID});
        FRET = T.FRETimageFA{rowID};
        %BW = logical(T.FRETimageFA{rowID});

    elseif hasFRETindex
        % Build binary mask from FRET index
        M = T.FRETindex{rowID};
        M(M == 100) = 0;         % set all bg 100 pixels to NaN
        FRET = M;
    else
    error('Row %d: neither FRETimageFA nor FRETindex found or non-empty.', rowID);
    end

    load('BlackJet.mat');
    nexttile;
    %FRET = T.FRETimageFA{rowID};
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
    
    
    % -------- Panel 3: T.ThicknessMap --------
    nexttile;

    %% Overlay image
    % TMap = T.ThicknessMap{rowID};   % 512×512 double, µm, 0 = background
    % 
    % if isempty(TMap)
    %     title('ThicknessMap (empty)');
    % else
    %     % Mask background (0) so it renders as white/blank
    %     TMapDisp = TMap;
    %     TMapDisp(TMap == 0) = NaN;         % NaN → not rendered by imagesc
    % 
    %     imagesc(TMapDisp);
    %     axis image off;
    %     colormap(gca, jet);                % thin=blue → thick=red
    %     cb = colorbar;
    %     cb.Label.String = 'Thickness (µm)';
    % 
    %     % CLim: 0 to 95th-percentile of foreground to avoid outlier stretch
    %     vals = TMap(TMap > 0);
    %     if ~isempty(vals)
    %         clim_max = prctile(vals, 95);
    %         if clim_max <= 0; clim_max = 1; end
    %         clim([0 clim_max]);
    %     end
    % 
    %     title('Junction Thickness (ThicknessMap)');
    % end

    %% Skeleton image 

    TMap = T.ThicknessMap{rowID};
    BW = TMap > 0;                  % foreground = any nonzero thickness
    Skel = bwskel(logical(BW));

    SkelDisp = TMap;
    SkelDisp(~Skel) = NaN;          % keep only skeleton pixels, NaN = not rendered

    imagesc(SkelDisp);
    axis image off;
    colormap(gca, hot);
    cb = colorbar;
    cb.Label.String = 'Thickness (µm)';

    vals = TMap(BW);
    if ~isempty(vals)
        clim_max = prctile(vals, 95);
        if clim_max <= 0; clim_max = 1; end
        caxis([0 clim_max]);
    end

    title('Junction Skeleton (ThicknessMap)');

end