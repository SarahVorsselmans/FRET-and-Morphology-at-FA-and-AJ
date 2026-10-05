function T = Calc_Thickness(T, Input)
% computeThickness: Calculates junction thickness and FRET/Intensity at skeleton
%
% Outputs added to T:
%   ThicknessMapSkel      - 512x512 double: thickness (µm) propagated to all fg pixels via skeleton
%   Thickness_Area        - cell: thickness (µm) at every foreground pixel
%   Thickness_Skeleton    - cell: thickness (µm) at every skeleton pixel
%   Thickness_Area_Av     - double: mean of Thickness_Area per row
%   Thickness_Skeleton_Av - double: mean of Thickness_Skeleton per row
%   Skel_FRET             - cell: mean FRET per skeleton pixel (via accumarray)
%   Skel_Int              - cell: mean intensity per skeleton pixel (via accumarray)
%   FRETavFiltered        - double: mean of Skel_FRET per row
%   Av_AverageIntensity   - double: mean of Skel_Int per row
%   skel_coords           - cell: [row, col] of skeleton pixels (for plotting)

    %% Parameters
    pxSize = Input.PixSize;   % µm per pixel
    pad    = 20;              % padding for EDT edge correction

    useFRET = ismember("FRETimageFA", T.Properties.VariableNames);
    useInt  = ismember("Int",         T.Properties.VariableNames);

    nRows = height(T);

    % Save memory
    Mask          = T.Mask;
    FRETimageFA   = T.FRETimageFA;
    Int           = T.Int;

    %% Pre-allocate all output columns
    T.FRETav                = nan(nRows,1);
    T.Av_AverageIntensity   = nan(nRows,1);
    T.ThicknessMapSkel      = cell(nRows, 1);
    % T.Thickness_Area      = cell(nRows, 1);
    T.Thickness_Skeleton    = cell(nRows, 1);
    T.Thickness_Area_Av     = nan(nRows, 1);
    T.Thickness_Skeleton_Av = nan(nRows, 1);
    T.Skel_FRET             = cell(nRows, 1);
    T.Skel_Int              = cell(nRows, 1);
    T.FRETavFiltered        = nan(nRows, 1);
    % T.FRETavSkel          = nan(nRows, 1);
    T.Av_AverageIntensityFiltered  = nan(nRows, 1);
    % T.Av_AverageIntensitySkel    = nan(nRows, 1);
    T.skel_coords           = cell(nRows, 1);

    fprintf('computeThickness: processing %d rows...\n', nRows);

    % Preallocate output arrays before the loop
    FRETav                         = nan(nRows, 1);
    FRETavFiltered                 = nan(nRows, 1);
    Av_AverageIntensity            = nan(nRows, 1);
    Av_AverageIntensityFiltered    = nan(nRows, 1);
    Thickness_Area_Av              = nan(nRows, 1);
    Thickness_Skeleton_Av          = nan(nRows, 1);
    % Cell arrays for variable-length outputs
    Skel_FRET                      = cell(nRows, 1);
    Skel_Int                       = cell(nRows, 1);
    Thickness_Skeleton             = cell(nRows, 1);
    ThicknessMapSkel               = cell(nRows, 1);

    %% Process each row
    for i = 1:nRows

        % Print progress every 50 rows so long runs stay observable.
        if mod(i, 50) == 0
            fprintf('  row %d / %d\n', i, nRows);
        end

        %% --- Raw image averages (no mask filtering applied) ---
        % Compute mean FRET and intensity directly from the raw images,
        % excluding NaN and 0, before any mask preprocessing is applied.
        if useFRET && ~isempty(FRETimageFA{i})
            rawFRET           = double(FRETimageFA{i}(:));         % all pixels as vector
            FRETav(i)         = mean(rawFRET(~isnan(rawFRET) & rawFRET ~= 0));
        end
        if useInt && ~isempty(Int{i})
            rawInt                      = double(Int{i}(:));         % all pixels as vector
            Av_AverageIntensity(i)      = mean(rawInt(~isnan(rawInt) & rawInt ~= 0));
        end

        %% --- Mask preprocessing ---
        BW = logical(Mask{i});
        BW = bwareaopen(BW, 100); % remove debris: discard connected components < 100 px
        BW = imclose(BW, strel('disk', 8)); % close small gaps / smooth contour with an 8-px disk

        if ~any(BW(:)) % guard: skip empty masks
            warning('Row %d: mask is empty after preprocessing, skipping.', i);
            continue
        end

        %% --- Padded Euclidean Distance Transform (EDT) ---
        % bwdist() assigns each pixel its distance to the nearest background pixel.
        % Padding prevents border pixels from being artificially close to the image
        % edge (which would underestimate thickness near the mask boundary).
        BWp = padarray(BW, [pad pad], 0, 'both');   % surround with 'pad' rows/cols of zeros
        Dp  = bwdist(~BWp);                         % EDT on padded mask: D(px) = dist to bg
        D   = Dp(1+pad : end-pad, 1+pad : end-pad); % strip padding → back to original size
        % D(px) now equals the radius of the largest inscribed disk centred at px.

        %% --- Morphological Skeleton ---
        % The skeleton is the medial axis: a 1-px-wide ridge running through the
        % centre of the adhesion. Thickness is most naturally sampled here because
        % skeleton pixels sit at the local maxima of D (the inscribed-disk radii).
        BWskel  = bwskel(BW);         % compute medial-axis skeleton of the binary mask
        skelIdx = find(BWskel);       % linear indices of all skeleton pixels (used throughout)
        nSkel   = numel(skelIdx);     % number of skeleton pixels

        if nSkel == 0 % guard: skip empty masks
            warning('Row %d: skeleton is empty, skipping.', i);
            continue
        end

        %% --- Nearest-skeleton mapping for ALL foreground pixels ---
        % For every foreground pixel we record which skeleton pixel is closest.
        % This lets us "propagate" skeleton-derived thickness out to the full mask.
        [~, idxMap] = bwdist(BWskel);  % idxMap(j) = linear index of the skeleton px nearest to px j
                                       % (second output of bwdist on a binary image)

        % Assign each skeleton pixel a unique integer label 1..nSkel.
        % skelLabel is zero everywhere except on the skeleton itself.
        skelLabel          = zeros(size(BW));
        skelLabel(skelIdx) = 1:nSkel;

        %% --- Compute thickness values ---
        % Local thickness = 2 × (inscribed-disk radius) = 2 × D, converted to µm.
        % Two variants are stored:
        %   • Skeleton version  – one value per skeleton pixel (direct, no propagation)
        %   • Area version      – one value per foreground pixel (each fg pixel inherits
        %                          the thickness of its nearest skeleton pixel)
    
        fgIdx   = find(BW);                         % linear indices of all foreground pixels
        groupID = skelLabel(idxMap(fgIdx));         % skeleton-pixel owner of each fg pixel
                                                    % (groupID used if per-branch stats are needed)
    
        % -- Area-weighted average thickness --
        % Each foreground pixel is assigned the diameter at its nearest skeleton pixel.
        tvals_fg_um        = 2 * D(idxMap(fgIdx)) * pxSize;   % propagated diameter [µm], one per fg px
        % T.Thickness_Area{i}  = tvals_fg_um;                 % (disabled) store full distribution
        Thickness_Area_Av(i) = mean(tvals_fg_um);           % scalar mean over all fg pixels
    
        % -- Skeleton-sampled thickness --
        % Thickness is read directly from D at each skeleton pixel (no propagation).
        % This avoids double-counting thick regions and is the primary output.
        tvals_skel_um              = 2 * D(skelIdx) * pxSize; % diameter [µm], one per skeleton px
        Thickness_Skeleton{i}    = tvals_skel_um;           % full distribution along skeleton
        Thickness_Skeleton_Av(i) = mean(tvals_skel_um);     % scalar mean along skeleton
    
        % -- Skeleton thickness map (512×512 spatial image) --
        % A zero image with thickness values painted onto skeleton pixels only.
        % Useful for visualising where thick/thin regions lie within the adhesion.
        thicknessMap          = zeros(size(BW));
        thicknessMap(skelIdx) = tvals_skel_um;       % only skeleton pixels are non-zero
        ThicknessMapSkel{i} = thicknessMap;
    
        % -- (Disabled) Area thickness map --
        % Alternative map where every foreground pixel carries its propagated thickness.
        % Re-enable if a filled (non-skeleton) spatial map is needed.
        % thicknessMap = zeros(size(BW));
        % thicknessMap(fgIdx) = tvals_fg_um;
        % T.ThicknessMapSkel{i} = thicknessMap;

        %% --- Skeleton pixel coordinates (for plotting) ---
        [row, col] = ind2sub(size(BW), skelIdx);
        T.skel_coords{i} = [row, col];

        %% --- FRET at skeleton (accumarray: mean FRET of all fg pixels per skel px) ---
        % For each skeleton pixel, compute the mean FRET value of all foreground
        % pixels that belong to it (i.e. whose nearest skeleton pixel is this one).
        % This gives a spatially resolved FRET map along the adhesion centreline.
        if useFRET && ~isempty(FRETimageFA{i})
            fgFRET         = double(FRETimageFA{i}(fgIdx)); % FRET values at all foreground pixels

            % accumarray groups fg pixels by their skeleton-pixel owner (groupID)
            % and computes the mean FRET within each group → one value per skeleton pixel.
            skel_FRET      = accumarray(groupID, fgFRET, [nSkel 1], @mean);
            Skel_FRET{i} = skel_FRET; % N×1 double: mean FRET per skeleton pixel

            % (Disabled) Alternative: average over skeleton pixels (propagated values).
            % FRETavSkel(i) = mean(skel_FRET);

            % Average over raw foreground pixels, excluding NaN and 0
            validFRET   = fgFRET(~isnan(fgFRET) & fgFRET ~= 0);
            FRETavFiltered(i) = mean(validFRET); % scalar mean FRET across valid fg pixels
        else
            warning('Row %d: no FRET data, skipping Skel_FRET.', i);
        end

        %% --- Intensity at skeleton (same logic as FRET) ---
        % Identical logic to the FRET block above, applied to the intensity image.
        % Produces a skeleton-resolved intensity map and a single scalar average.
        if useInt && ~isempty(Int{i})
            fgInt         = double(Int{i}(fgIdx));
            skel_Int      = accumarray(groupID, fgInt, [nSkel 1], @mean);
            Skel_Int{i} = skel_Int;

            % (Disabled) Alternative: average over skeleton pixels (propagated values).
            % Av_AverageIntensitySkel(i) = mean(skel_Int);

            % Average from raw fg pixels, excluding NaN and 0
            validInt                         = fgInt(~isnan(fgInt) & fgInt ~= 0);
            Av_AverageIntensityFiltered(i) = mean(validInt);
        else
            warning('Row %d: no Int data, skipping Skel_Int.', i);
        end
    end
    % Write all results to T in one go after the loop
    T.FRETav                  = FRETav;
    T.Av_AverageIntensity     = Av_AverageIntensity;
    T.FRETavFiltered                  = FRETavFiltered;
    T.Av_AverageIntensityFiltered     = Av_AverageIntensityFiltered;
    T.Thickness_Area_Av       = Thickness_Area_Av;
    T.Thickness_Skeleton_Av   = Thickness_Skeleton_Av;
    T.Skel_FRET               = Skel_FRET;
    T.Skel_Int                = Skel_Int;
    T.Thickness_Skeleton      = Thickness_Skeleton;
    T.ThicknessMapSkel        = ThicknessMapSkel;

    fprintf('Done.\n');
end


% function T = computeThickness(T, Version)
% % computeThickness: Calculates junction thickness and plots distribution
% % Inputs:
% %   T.Mask    - logical 2D array (segmentation mask)
% %   Version   - 'Skeleton' or 'Area'
% % Outputs:
% %   T.ThicknessAv  - mean thickness (µm)
% %   T.ThicknessSum - table with summary stats + all values + Version label
% 
%     useFRET = ismember("FRETimageFA", T.Properties.VariableNames);
% 
%     %% Parameters
%     pxSize      = 0.24;    % µm per pixel (adjust if needed)
%     pad         = 20;      % padding for EDT (px)
%     %binWidth_um = 0.5; % bin width for histogram-line plot (µm)
% 
%     %% Prepare output columns
%     nRows = height(T);
%     T.ThicknessMap = cell(nRows, 1);
% 
%     fieldAv  = ['ThicknessAv_' Version];
%     %fieldSum = ['ThicknessSum_' Version];
%     fieldAll = ['Thickness_' Version];
% 
%     T.(fieldAv)  = nan(nRows,1);
%     %T.(fieldSum) = cell(nRows,1);  % each cell will contain a 1x1 table
%     T.(fieldAll) = cell(nRows,1);
% 
%     %% Process each row
%     for i = 1:nRows
%         %% Choose source of mask
%         BW = T.Mask{i};
% 
%         %% Preprocessing
%         %BW = logical(T.Mask);
%         %BW = T.Mask;
%         BW = bwareaopen(BW, 100);     % remove small objects (< 100 px)
%         % BW = imfill(BW, 4, 'holes');  % fill holes (not diagonal (4))
%         % BW = imclose(BW, strel('disk', 1)); % close tiny gaps
%         % BW = bwareaopen(~BW, 500, 4); % remove small holes (< 500 px)
%         % BW = ~BW;                     % invert back
% 
%         %BW_original = BW;
%         % Use closing to find regions that SHOULD be filled
%         SE = strel('disk', 8);
%         BW = imclose(BW, SE);
%         %imagesc(BW_closed)
%         % Find holes in the closed version
%         % holes = ~BW;
%         % holes = bwareaopen(holes, 500, 4); % keep large true background
%         % BW_closed_filled = ~holes;  % closed version with small holes filled
%         % imagesc(BW_closed_filled)
%         % Apply ONLY those filled regions to the original
%         % BW = BW_original | BW_closed_filled;
%         %imagesc(BW)
% 
%         %% Skeleton 
%         BWskel = bwskel(BW);
% 
%         %% Euclidean Distance Transform with PADDING
%         BWp = padarray(BW, [pad pad], 0, 'both'); % pad with background (false)
%         Dp  = bwdist(~BWp); % EDT on padded image
%         D   = Dp(1+pad:end-pad, 1+pad:end-pad); % crop back
% 
%         %% Thickness calculation
%         thickness_area = 2 * D; % radius -> diameter in pixels
%         thickness_skeleton = thickness_area .* BWskel;
% 
%         % Propagate thickness to all foreground pixels
%         % Why?: That way the value of every pixel in e.g. a thick section has a thick value and not only the pixels at the center :)
%         [~, idxMap] = bwdist(BWskel); % nearest skeleton pixel index
%         thickness_map = zeros(size(BW));
%         thickness_map(BW) = thickness_area(idxMap(BW)); % area-weighted field (px)
% 
%         thickness_map_um = thickness_map * pxSize;
%         thickness_map_um(~BW) = 0;          % background stays 0
%         T.ThicknessMap{i} = thickness_map_um;
% 
%         %% Select version
%         switch Version
%             case 'Skeleton' % — each position along the junction counts once (thinner weigh higher?)
%                 tvals_px = thickness_skeleton(BWskel);
%                 tvals_px = tvals_px(tvals_px > 0);
%             case 'Area' % — thicker regions dominate
%                 tvals_px = thickness_map(BW);
%             otherwise
%                 error('Version must be ''Skeleton'' or ''Area''.');
%         end
% 
%         %% Convert to µm
%         tvals_um = tvals_px * pxSize;
% 
%         %% Summary stats
%         meanVal = mean(tvals_um);
%         % medianVal = median(tvals_um);
%         % stdVal = std(tvals_um);
%         % minVal = min(tvals_um);
%         % maxVal = max(tvals_um);
%         % nVal = numel(tvals_um);
% 
%         % Save outputs
%         T.(fieldAv)(i) = meanVal;
%         T.(fieldAll){i} = tvals_um;
% 
%         % summaryStruct = struct( ...
%         %     'Mean_um',     meanVal, ...
%         %     'Median_um',   medianVal, ...
%         %     'Std_um',      stdVal, ...
%         %     'Min_um',      minVal, ...
%         %     'Max_um',      maxVal, ...
%         %     'N',           nVal, ...
%         %     'AllValues_um', {tvals_um(:)}, ...     % <-- wrap in a cell to keep it 1x1
%         %     'Version',     string(Version) );
%         % T.(fieldSum){i} = struct2table(summaryStruct, 'AsArray', true);
% 
%         %% Plot histogram as line plot (bin size = 0.5 µm)
%         % Bin edges common and neat
%         % tmax       = ceil(max(tvals_um)/binWidth_um) * binWidth_um;
%         % edges_um   = 0:binWidth_um:tmax;
%         % centers_um = edges_um(1:end-1) + diff(edges_um)/2;
%         % 
%         % % Fractions per bin = PDF * bin_width
%         % fractions = histcounts(tvals_um, edges_um, 'Normalization','pdf') .* diff(edges_um);
% 
%         % figure('Name', sprintf('Thickness Distribution — Row %d (%s)', i, Version), 'Color','w');
%         % hold on;
%         % plot(centers_um, fractions, '-o', 'LineWidth', 2, 'MarkerSize', 4, 'Color', [0 0.45 0.74]);
%         % fill([centers_um fliplr(centers_um)], [fractions zeros(size(fractions))], ...
%         %     [0 0.45 0.74], 'FaceAlpha', 0.2, 'EdgeColor','none');
%         % xlabel('Junction thickness (\mum)');
%         % ylabel('Fraction of pixels');
%         % title(sprintf('Thickness distribution — Row %d (%s)', i, Version));
%         % grid on; box on;
%         % hold off;
%     end
% end
