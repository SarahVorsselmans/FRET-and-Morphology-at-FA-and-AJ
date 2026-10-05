function T = Calc_FA_AllMetrics(T, Input, options)
% Inputs:
%   T      - Table with columns:
%               .Int{i}   - intensity image (double)
%               .Mask{i}  - binary mask (preferred source)
%             Optional fallbacks: .FRETimageFA{i}, .FRETindex{i}
%   Input  - Struct with fields:
%               .PixSize  - pixel size in µm (e.g. 0.16)
%
% Output columns added to T:
%   .CC_Obj
%   .Area, .Ecc, .MajAx, .MinAx, .MinMajRatio, .Orientations
%   .Centroids, .NND, .NND3, .FADens_r5
%   .IntensityPerPixel, .TotalInt, .AverageIntensity, .AverageIntensityNorm, .MaxInt, .MinInt

% ── Inputs ────────────────────────────────────────
arguments
    T                table
    Input            struct
    options.FRET     (1,1) logical = false
end
FRET = options.FRET;

% ── Detect available mask columns ────────────────────────────────────────────
    hasMask = ismember("Mask", T.Properties.VariableNames);
    hasInt = ismember("Int",   T.Properties.VariableNames);
    if FRET
        hasFRET = ismember("FRETimageFA", T.Properties.VariableNames);
    else
        hasFRET = 0;
    end

    if ~hasMask || ~hasInt || ~hasFRET
        error('Calc_FA_AllMetrics: No mask, int (or FRET) source found.');
    end

    PixSize   = Input.PixSize;
    radius_um = 5;          % fixed radius for FA density (µm)
    nRows     = height(T);

% ── Preallocate all output columns ───────────────────────────────────────────
    % All [Nx1 double]!
    % Basic metrics
    T.CC_Obj       = cell(nRows,1);
    if FRET
        T.FRET     = cell(nRows,1);
        T.FRETav   = nan(nRows,1);
    end
    T.Av_AverageIntensity = nan(nRows,1);

    % Shape metrics
    T.Area         = cell(nRows,1); % - Area per FA in pixels
    T.Ecc          = cell(nRows,1); % - Eccentricity per FA ('how stretched the ellipse is')
    T.MajAx        = cell(nRows,1); % - Long axis of fitted ellipse per FA (length)
    T.MinAx        = cell(nRows,1); % - Short axis of fitted ellipse per FA (width)
    T.MinMajRatio  = cell(nRows,1); % - Roundness per FA (1 = round | ~0 = line)
    T.Orientations = cell(nRows,1); % - Orientation per FA (-90 to +90 degrees)

    % Clustering metrics
    T.Centroids  = cell(nRows,1); % - Nx2 centroid coordinates (pixels) (but cents_um is in µm)
    T.NND        = cell(nRows,1); % - Nearest-neighbor distance per FA (µm)
    T.NND3       = cell(nRows,1); % - 3rd nearest neighbor distance per FA (µm)
    T.FADens_r5  = cell(nRows,1); % - FA density within 5 µm radius per FA (FA/µm²)

    % Per-pixel intensity metrics
    T.IntensityPerPixel    = cell(nRows,1);   % - Average intensity per pixel per FA
    T.TotalInt             = cell(nRows,1);   % - Summed intensity per FA
    T.AverageIntensity     = cell(nRows,1);   % - Mean intensity per FA
    T.AverageIntensityNorm = cell(nRows,1);   % - Z-score normalized mean int
    T.MaxInt               = cell(nRows,1);   % - Max intensity per FA
    T.MinInt               = cell(nRows,1);   % - Min intensity per FA

    % Save memory
    Mask          = T.Mask;
    FRETimageFA   = T.FRETimageFA;
    Int           = T.Int;

% ── Main loop ────────────────────────────────────────────────────────────────
    for i = 1:nRows

        % ── 1. Resolve binary mask ──
        BW = Mask{i}; % already binary

        % ── 2. Connected components ─────────────────────────────
        CC    = bwconncomp(BW, 8);
        nFA   = CC.NumObjects;

        % Store CC object (used downstream if needed)
        T.CC_Obj{i} = CC;

        % ── 3. Handle empty images ────────────────────────────────────────────
        if nFA == 0
            T.Area{i}        = [];
            T.Ecc{i}         = [];
            T.MajAx{i}       = [];
            T.MinAx{i}       = [];
            T.MinMajRatio{i} = [];
            T.Orientations{i}= [];
            T.Centroids{i}   = [];
            T.NND{i}         = [];
            T.NND3{i}        = [];
            T.FADens_r5{i}   = [];
            T.IntensityPerPixel{i}    = [];
            T.TotalInt{i}             = [];
            T.AverageIntensity{i}     = [];
            T.AverageIntensityNorm{i} = [];
            T.MaxInt{i}      = [];
            T.MinInt{i}      = [];
            if FRET
                T.FRET{i}    = [];
                T.FRETav(i)  = NaN;
            end
            T.Av_AverageIntensity(i) = NaN;
            continue
        end

        % ── 4. regionprops ──────────────────────────────────────
        img                 = double(Int{i});
        validInt            = img(~isnan(img) & img ~= 0);
        Av_AverageIntensity = mean(validInt);
        if FRET
            fret        = double(FRETimageFA{i});
            % extra:calc FRETav as usual
            validFRET   = fret(~isnan(fret) & fret ~= 0);
            FRETav      = mean(validFRET);
        end
        stats = regionprops(CC, ...
            'Area','Eccentricity', ...
            'MajorAxisLength','MinorAxisLength', ...
            'Orientation','Centroid');

        % ── 5. Per-FA loop ──────────────────────────────────────
        totInt   = zeros(nFA,1);
        avgInt   = zeros(nFA,1);
        maxInt   = zeros(nFA,1);
        minInt   = zeros(nFA,1);
        intPerPx = zeros(nFA,1);
        cents_px = nan(nFA,2);
        if FRET
            avgfret = zeros(nFA,1);
        end

        for j = 1:nFA
            pix  = CC.PixelIdxList{j};

            % Int
            vals = img(pix);
            totInt(j)   = sum(vals);
            avgInt(j)   = mean(vals);
            maxInt(j)   = max(vals);
            minInt(j)   = min(vals);
            intPerPx(j) = totInt(j) / numel(pix);
            
            % FRET
            if FRET
                valsfret   = fret(pix);
                avgfret(j) = mean(valsfret);
            end

            % Centroid (from regionprops, no recompute needed)
            if ~isempty(stats(j).Centroid)
                cents_px(j,:) = stats(j).Centroid;
            end
        end

        % ── 6. Shape metrics (vectorised from regionprops) ───────────────────
        Ecc   = vertcat(stats.Eccentricity);
        Area  = vertcat(stats.Area);
        MajAx = vertcat(stats.MajorAxisLength);
        MinAx = vertcat(stats.MinorAxisLength);
        MinMajRatio = MinAx ./ MajAx;

        % Orientation — suppress unreliable estimates for near-circular FAs
        aspectRatio = MajAx ./ max(MinAx, 1);
        keep        = aspectRatio > 1.5; % Set to 2 for min/maj ratio of 0.5
        angles_deg  = nan(nFA,1);
        rawAngles   = vertcat(stats.Orientation);
        angles_deg(keep) = rawAngles(keep);

        % ── 7. Clustering metrics ─────────────────────────────────────────────
        % 1) Nearest Neighbor Distance (NND)
        % Distances between centroids
        cents_um = cents_px * PixSize;
        D  = squareform(pdist(cents_um));               % NxN pairwise distances (µm)
        N  = size(D,1);
        % Removing self-distances: diagonal entries are zero (distance to itself)
        D2 = D;
        D2(1:N+1:end) = inf;                            % Setting them to Inf ensures they won’t be picked as nearest neighbors
        % Sorting distances per object
        sortedD = sort(D2, 2, 'ascend');                % -> sortedD(i,1) = 1st neigbor - sortedD(i,2) = 2nd nearest neigbor etc.
        D(D == Inf) = NaN; % or just exclude those pairs

        if size(sortedD, 1) < 2 || size(sortedD, 2) < 2
            %warning('Category too small for NND calculation, skipping.');
            NND = NaN;
        else        
            NND  = sortedD(:,min(1,end));               % guard if nFA < 1
            NND3 = sortedD(:,min(3,end));               % guard if nFA < 3
        end

        % 2) Amount of neigbors in radius 
        % Find neighbors within a radius
        withinR    = (D <= radius_um) & (D > 0);        % Logical matrix: true = object j is within radius of object i & excludes self (D > 0)
        neighCount = sum(withinR, 2);                   % Count neighbors
        FADens_r5  = neighCount / (pi * radius_um^2);   % Convert to density (obj/µm²)

        % ── 8. Intensity metrics ──────────────────────────────────────────────
        mu    = mean(avgInt);
        sigma = std(avgInt);
        if sigma > 0
            avgIntNorm = (avgInt - mu) / sigma;
        else
            avgIntNorm = zeros(nFA,1);
        end

        % ── 9. Write results to table ─────────────────────────────────────────
        % Fret
        if FRET
            T.FRET{i}    = avgfret;
            T.FRETav(i)  = FRETav;
        end
        T.Av_AverageIntensity(i) = Av_AverageIntensity;
        
        % Shape
        T.Area{i}        = Area;
        T.Ecc{i}         = Ecc;
        T.MajAx{i}       = MajAx;
        T.MinAx{i}       = MinAx;
        T.MinMajRatio{i} = MinMajRatio;
        T.Orientations{i}= angles_deg;

        % Clustering
        T.Centroids{i}  = cents_px;
        T.NND{i}        = NND;
        T.NND3{i}       = NND3;
        T.FADens_r5{i}  = FADens_r5;

        % Intensity
        T.IntensityPerPixel{i}    = intPerPx;
        T.TotalInt{i}             = totInt;
        T.AverageIntensity{i}     = avgInt;
        T.AverageIntensityNorm{i} = avgIntNorm;
        T.MaxInt{i}               = maxInt;
        T.MinInt{i}               = minInt;

    end % row loop
end