function T = Calc_FA_AveragingMetrics(T, Input, options)
% Calc_FA_Metrics
% -------------------------------------------------------------
% Computes FA metrics in both pixel and micron units.
%
% Computes:
%   • FRETav2               (mean FRET per FA (if FRET = true))
%   • FA_Count              (# of adhesion sites)
%   • Av_PixSize            (mean FA size in pixels)
%   • Total_FA_Area_pix     (sum of FA pixel areas)
%   • Median_FA_Size_pix    (median FA pixel area)
%   • Av_UmSize             (mean FA size in µm²)
%   • Total_FA_Area_um      (sum of FA areas in µm²)
%   • Av_MajAx_um           (mean major axis length in µm)
%   • Av_MinMajRatio        (mean minor/major axis ratio)
%   • Av_Ecc                (mean eccentricity)
%   • Av_NND3               (mean nearest-neighbour distance)
%   • Av_FADens_r5          (mean FA density, radius 5)
%   • Av_MajAx              (mean major axis in pixels)
%   • Av_AverageIntensity2   (mean per-FA average intensity)
%   • MeanAngle             (mean resultant angle)
%   • CircVar               (circular variance [0=aligned, 1=random])
%   • AlignIndex            (1 - CircVar, easier to read)
%
% INPUT:
%   T          - table with cell-array columns: Area, MinMajRatio, Ecc,
%                NND3, FADens_r5, MajAx, AverageIntensity
%   PixSize_um - scalar, pixel size in µm (used for unit conversion)
%                Area columns scale by PixSize_um^2; length columns by PixSize_um
%
% OUTPUT:
%   T with all new metric columns appended.
%   Only calculates mean if there are minimum of 5 objects
%
% Notes:
%   - Does NOT reorder existing columns.
%   - Safely handles: empty cells, non-numeric entries,
%     row/column vectors, and NaNs within vectors.
% -------------------------------------------------------------

% ── Inputs ────────────────────────────────────────
arguments
    T                table
    Input            struct
    options.FRET     (1,1) logical = false
end
FRET        = options.FRET;
PixSize_um  = Input.PixSize;

    %% --- Input validation -------------------------------------------
    if ~ismember('Area', T.Properties.VariableNames)
        warning('Calc_FA_Metrics:missingColumn', ...
            'Area column not found. Returning T unchanged.');
        return;
    end
    
    %% --- Column names to pull from T --------------------------------
    colNames = {'Area','MinMajRatio','Ecc','NND3','FADens_r5','MajAx','AverageIntensity'};
    
    % Warn once per missing optional column rather than erroring mid-loop
    missingCols = colNames(~ismember(colNames, T.Properties.VariableNames));
    if ~isempty(missingCols)
        warning('Calc_FA_Metrics:missingColumn', ...
            'Missing columns (will be NaN): %s', strjoin(missingCols, ', '));
    end
    
    %% --- Preallocate ------------------------------------------------
    n = height(T);
    
    if FRET
        FRETav2         = nan(n,1);
    end
    FA_Count            = nan(n,1);
    Av_PixSize          = nan(n,1);
    Total_FA_Area_pix   = nan(n,1);
    Median_FA_Size_pix  = nan(n,1);
    Av_UmSize           = nan(n,1);   % µm
    Total_FA_Area_um    = nan(n,1);   % µm
    Av_MinMajRatio      = nan(n,1);
    Av_Ecc              = nan(n,1);
    Av_NND3             = nan(n,1);   % µm
    Av_FADens_r5        = nan(n,1);   % µm
    Av_MajAx            = nan(n,1);
    Av_MajAx_um         = nan(n,1);   % µm
    Av_AverageIntensity2 = nan(n,1);
    Av_MeanAngle        = nan(n,1);
    Av_CircVar          = nan(n,1);
    Av_AlignIndex       = nan(n,1);
    
    %% --- Helper: safely extract a numeric column-vector from a cell --
        function v = getVec(tbl, col, row)
            if ~ismember(col, tbl.Properties.VariableNames)
                v = []; return;
            end
            raw = tbl.(col){row};
            if isnumeric(raw) && isvector(raw) && ~isempty(raw)
                v = raw(:);          % guarantee column vector
            else
                v = [];
            end
        end
    
    %% --- Main loop --------------------------------------------------
    for i = 1:n
    
        v  = getVec(T, 'Area',            i);  if isempty(v), continue; end
        v2 = getVec(T, 'MinMajRatio',     i);
        v3 = getVec(T, 'Ecc',             i);
        v4 = getVec(T, 'NND3',            i);
        v5 = getVec(T, 'FADens_r5',       i);
        v6 = getVec(T, 'MajAx',           i);
        v7 = getVec(T, 'AverageIntensity',i);
        v8 = getVec(T, 'Orientations',    i);
        if FRET
            v9 = getVec(T, 'FRET', i);
        end

        % --- Pixel-space metrics ---
        FA_Count(i)           = sum(~isnan(v)); % count non-NaN FAs

        if FA_Count(i) >= 5
            Av_PixSize(i)         = mean(v,  'omitnan');
            Total_FA_Area_pix(i)  = sum(v,   'omitnan');
            Median_FA_Size_pix(i) = median(v,'omitnan');
            
            % --- FRET metrics ---
            if FRET
                FRETav2(i)        = mean(v9,  'omitnan');
            end
    
            % --- µm conversions (area: px² → µm²;  length: px → µm) ---
            Av_UmSize(i)          = Av_PixSize(i)        * PixSize_um^2;
            Total_FA_Area_um(i)   = Total_FA_Area_pix(i) * PixSize_um^2;
    
            % --- Orientation metrics ---
            valid = ~isnan(v8);
            % Convert to radians, double to map 180° -> 360° circle
            theta = deg2rad(v8(valid)) * 2;
            % Mean resultant vector
            C = mean(cos(theta));
            S = mean(sin(theta));
            R = sqrt(C^2 + S^2); % mean resultant length [0,1]
            % Unwrap mean angle back to [0, 180)
            meanAngle = rad2deg(atan2(S, C)) / 2;
            if meanAngle < 0
                meanAngle = meanAngle + 180;
            end
            Av_MeanAngle(i) = meanAngle;
            Av_CircVar(i)   = 1 - R;
            Av_AlignIndex(i)= R;
    
            % --- Shape / spatial metrics ---
            if ~isempty(v2), Av_MinMajRatio(i)      = mean(v2,'omitnan'); end
            if ~isempty(v3), Av_Ecc(i)              = mean(v3,'omitnan'); end
            if ~isempty(v4), Av_NND3(i)             = mean(v4,'omitnan'); end
            if ~isempty(v5), Av_FADens_r5(i)        = mean(v5,'omitnan'); end
            if ~isempty(v6)
                Av_MajAx(i)                         = mean(v6,'omitnan');
                Av_MajAx_um(i)                      = Av_MajAx(i) * PixSize_um;
            end
            if ~isempty(v7), Av_AverageIntensity2(i) = mean(v7,'omitnan'); end

        else
            if FRET
                FRETav2(i)         = NaN;
            end
            Av_PixSize(i)          = NaN;
            Av_UmSize(i)           = NaN;
            Median_FA_Size_pix(i)  = NaN;
            Total_FA_Area_pix(i)   = NaN;
            Total_FA_Area_um(i)    = NaN;
            Av_MinMajRatio(i)      = NaN;
            Av_Ecc(i)              = NaN;
            Av_NND3(i)             = NaN;
            Av_FADens_r5(i)        = NaN;
            Av_MajAx(i)            = NaN;
            Av_MajAx_um(i)         = NaN; 
            Av_AverageIntensity2(i) = NaN;
            Av_MeanAngle(i)        = NaN;
            Av_CircVar(i)          = NaN;
            Av_AlignIndex(i)       = NaN;
        end
    end
    
    %% --- Attach to table --------------------------------------------
    if FRET
        T.FRETav2         = FRETav2;
    end
    T.FA_Count            = FA_Count;
    T.Av_PixSize          = Av_PixSize;
    T.Av_UmSize           = Av_UmSize;
    T.Median_FA_Size_pix  = Median_FA_Size_pix;
    T.Total_FA_Area_pix   = Total_FA_Area_pix;
    T.Total_FA_Area_um    = Total_FA_Area_um;
    T.Av_MinMajRatio      = Av_MinMajRatio;
    T.Av_Ecc              = Av_Ecc;
    T.Av_NND3             = Av_NND3;
    T.Av_FADens_r5        = Av_FADens_r5;
    T.Av_MajAx            = Av_MajAx;
    T.Av_MajAx_um         = Av_MajAx_um; 
    T.Av_AverageIntensity2 = Av_AverageIntensity2;
    T.Av_MeanAngle        = Av_MeanAngle;
    T.Av_CircVar          = Av_CircVar;
    T.Av_AlignIndex       = Av_AlignIndex;

end

% % =====================================================
% %% HELPER to saveguard averaging of only >= 5 objects
% % =====================================================
% 
% function m = safeMean(x, minN)
%     valid = ~isnan(x);
%     if sum(valid) >= minN
%         m = mean(x(valid));
%     else
%         m = NaN;
%     end
% end