function [mask] = ROIandSeg(data, flags, minSize)

% Input: data.mask = bg is black, seg is white
if ~isfield(flags, 'ContinueWithoutROI')
    flags.ContinueWithoutROI = 0;
end

disp(['--- Analyzing ' data.sampleIndex ' ---']) 

%% Extracting Data
if isfield(data, 'int')
    Int = data.int;
else
    Int = data.DD;
end

%% ROI Selection and Adjusted Mask
% Step 1: clean the original mask
if isfield(data, 'mask') % check if this variable exists
    SizeOrientation = 8; % 8: looks at surrounding pixels including diagonal; 4: not diagonal
    cleanedMask = bwareaopen(data.mask, minSize, SizeOrientation);  % small objects removed % Obj = 1, BG = 0
    smallRegions = data.mask & ~cleanedMask; % Obj = 1 (we don't want these Obj)
else
    % IllastikMask = []; % creates an empty matrix just to make sure the code doesn't error
end

if ismember(flags.manualROI, [1, 2])
    % Ensure ROI folder exists
    ROIpath = fullfile(data.dataPath, 'ROI');
    if ~exist(ROIpath, 'dir')
        mkdir(ROIpath);
    end

    % Determine ROI file name and inversion
    if flags.manualROI == 1
        ROIfile = fullfile(ROIpath, sprintf('ROIRemove_%s.mat', data.name));
        invertROI = false;
    else
        ROIfile = fullfile(ROIpath, sprintf('ROIKeep_%s.mat', data.name));
        invertROI = true;
    end

    % Flag - flags.ContinueWithoutROI
    skipROI = exist(ROIfile, 'file') == 0 && ...
              flags.ContinueWithoutROI == 1;
    
    if skipROI
        % --- FALLBACK TO "NO ROI" BEHAVIOR ---
        mask = cleanedMask;
        mask_adj = mask;
    
        % Save adjusted mask as TIFF (same as original else branch)
        tiffFilename = fullfile(data.dataPath, sprintf('Cell%d_AdjSeg.tiff', data.index));
        imwrite(uint8(mask_adj)*255, tiffFilename);
        fprintf('ContinueWithoutROI.\n');
        return; % ! exit function here
    end

    % Step 2: overlay for interactive drawing
    %overlayMask = ~cleanedMask; 
    overlayMask = cleanedMask; % objects=1, background=0

    % Step 3: load existing ROI or draw new
    if exist(ROIfile, 'file') && ~isfield(flags, 'additionalmanualROI') || ...
       exist(ROIfile, 'file') && isfield(flags, 'additionalmanualROI') && ~flags.additionalmanualROI
        % Normal behaviour: load and use existing ROI, skip interactive
        tmp = load(ROIfile);
        RegionMask = tmp.RegionMask;
    else
        % No existing ROI, OR additionalmanualROI=1: run interactive
        if exist(ROIfile, 'file') && isfield(flags, 'additionalmanualROI') && flags.additionalmanualROI
            % Load existing ROI and pre-apply it to the overlay so the user
            % sees what has already been removed
            tmp = load(ROIfile);
            existingROI = tmp.RegionMask;
            overlayMask = cleanedMask & ~existingROI; % show already-removed regions as gone
            overlayMask = bwareaopen(overlayMask, minSize, SizeOrientation);
            fprintf('Loaded existing ROI for %s — draw additional removals.\n', data.name);
        else
            existingROI = false(size(cleanedMask));
        end
    
        % Draw additional (or first) ROIs
        newROI = interactive_roi(Int, overlayMask);
    
        if invertROI
            newROI = ~newROI;
        end
    
        % Merge existing ROI + new small fragments + newly drawn ROI
        newSmallRegions = (cleanedMask & ~existingROI) & ~overlayMask;
        RegionMask = existingROI | newSmallRegions | newROI;
    end

    % % OLD - Step 3: load existing ROI or draw new
    % if exist(ROIfile, 'file')
    %     tmp = load(ROIfile);
    %     RegionMask = tmp.RegionMask; % Obj (to be removed) = 1, Bg = 0
    % else
    %     RegionMask = interactive_roi(Int, overlayMask);
    %     if invertROI % if true
    %         RegionMask = ~RegionMask;
    %     end
    % end

    % Step 4: combine ROI with small excluded regions
    RegionMask = RegionMask | smallRegions; % Obj (small & ROI) = 1

    % Step 5: create adjusted mask by removing ROI from original
    mask_adj = cleanedMask & ~RegionMask; % Obj(seg) = 1 

    % Step 7: save ROI as .mat
    save(ROIfile, 'RegionMask');

    % Step 8: save adjusted mask as TIFF
    tiffFilename = fullfile(ROIpath, sprintf('Cell%d_AdjSeg.tiff', data.index));
    imwrite(uint8(mask_adj)*255, tiffFilename);

    % Output
    mask = mask_adj;

else
    % No manual ROI: just clean original mask
    mask = cleanedMask;
    mask_adj = mask;

    % Save adjusted mask as TIFF
    tiffFilename = fullfile(data.dataPath, sprintf('Cell%d_AdjSeg.tiff', data.index));
    imwrite(uint8(mask_adj)*255, tiffFilename);
end