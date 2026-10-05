function T = Calc_FA_TensionPerDensity(T, options)
% Calc_FA_TensionPerDensity
% -------------------------------------------------------------
% Computes a per-FA "tension efficiency" ratio (tension / density)
% and aggregates it per row (image/cell), following the same
% cell-array-of-vectors structure used by Calc_FA_AllMetrics /
% Calc_FA_AveragingMetrics.
%
% Computes (per row):
%   • TensionPerDens        - Nx1 cell, per-FA ratio (tension/density)
%   • TensionPerDens_log10  - Nx1 cell, log10 of the per-FA ratio
%   • Av_TensionPerDens     - mean per-FA ratio (row average, >=5 FAs)
%   • Med_TensionPerDens    - median per-FA ratio (row average, >=5 FAs)
%   • Av_TensionPerDens_log10 - mean of log10 ratio (row average)
%   • N_TensionPerDens      - number of FAs that passed the density filter
%
% INPUT:
%   T                - table with cell-array columns for tension and density
%                      (default: 'FRET' and 'AverageIntensity', one value per FA)
%   options.TensionCol  - name of per-FA tension column   (default 'FRET')
%   options.DensityCol  - name of per-FA density column   (default 'AverageIntensity')
%   options.MinDensity  - minimum density value to include an FA in the
%                         ratio (filters out zero/near-zero denominators)
%                         (default: 0, i.e. only strictly positive values kept)
%   options.MinCount    - minimum number of valid FAs per row required to
%                         compute a row average (default 5, matches
%                         Calc_FA_AveragingMetrics convention)
%   options.ScaleFactor - multiplicative constant applied to the ratio
%                         purely for readability (e.g. 1000 to turn
%                         0.003 into 3.0). Does NOT change any statistics
%                         (means/SDs/t-tests scale identically), only the
%                         reported units. Default: 1 (no scaling).
%
% OUTPUT:
%   T with new columns appended (does not reorder existing columns).
%
% Notes:
%   - FAs with density <= options.MinDensity (or NaN tension/density) are
%     excluded from the ratio for that FA (set to NaN in the per-FA output).
%   - Row-level averages are only computed if >= options.MinCount valid
%     FAs remain, otherwise NaN (same convention as Calc_FA_AveragingMetrics).
%   - ScaleFactor is applied before the log10 transform, so
%     TensionPerDens_log10 = log10(ScaleFactor * ratio). Since
%     log10(k*x) = log10(k) + log10(x), scaling shifts the log10 values
%     by a constant (log10(ScaleFactor)) but does not change differences
%     between conditions, and therefore does not affect stats on the
%     log-transformed values either.
% -------------------------------------------------------------

arguments
    T                   table
    options.TensionCol  (1,1) string = "FRET"
    options.DensityCol  (1,1) string = "AverageIntensity"
    options.MinDensity  (1,1) double = 5
    options.MinCount    (1,1) double = 5
    options.ScaleFactor (1,1) double = 1
end

TensionCol  = options.TensionCol;
DensityCol  = options.DensityCol;
MinDensity  = options.MinDensity;
MinCount    = options.MinCount;
ScaleFactor = options.ScaleFactor;

%% --- Input validation -------------------------------------------
if ~ismember(TensionCol, T.Properties.VariableNames) || ...
   ~ismember(DensityCol, T.Properties.VariableNames)
    error('Calc_FA_TensionPerDensity:missingColumn', ...
        'Required columns "%s" and/or "%s" not found in T.', TensionCol, DensityCol);
end

%% --- Preallocate --------------------------------------------------
n = height(T);

T.TensionPerDens       = cell(n,1);
T.TensionPerDens_log10 = cell(n,1);
Av_TensionPerDens        = nan(n,1);
Med_TensionPerDens       = nan(n,1);
Av_TensionPerDens_log10  = nan(n,1);
N_TensionPerDens         = nan(n,1);

%% --- Helper: safely extract a numeric column-vector from a cell ---
    function v = getVec(tbl, col, row)
        raw = tbl.(col){row};
        if isnumeric(raw) && isvector(raw) && ~isempty(raw)
            v = raw(:);
        else
            v = [];
        end
    end

%% --- Main loop ------------------------------------------------------
for i = 1:n

    tension = getVec(T, TensionCol, i);
    density = getVec(T, DensityCol, i);

    if isempty(tension) || isempty(density) || numel(tension) ~= numel(density)
        T.TensionPerDens{i}       = [];
        T.TensionPerDens_log10{i} = [];
        continue
    end

    % ── Filter: valid = non-NaN tension, density above threshold ──
    valid = ~isnan(tension) & ~isnan(density) & (density > MinDensity);

    ratio = nan(size(tension));
    ratio(valid) = ScaleFactor * (tension(valid) ./ density(valid));

    ratio_log10 = nan(size(ratio));
    posRatio = valid & ratio > 0;
    ratio_log10(posRatio) = log10(ratio(posRatio));

    T.TensionPerDens{i}       = ratio;
    T.TensionPerDens_log10{i} = ratio_log10;

    % ── Row-level aggregation (only if enough valid FAs) ──
    nValid = sum(valid);
    N_TensionPerDens(i) = nValid;

    if nValid >= MinCount
        Av_TensionPerDens(i)       = mean(ratio(valid), 'omitnan');
        Med_TensionPerDens(i)      = median(ratio(valid), 'omitnan');
        Av_TensionPerDens_log10(i) = mean(ratio_log10(posRatio), 'omitnan');
    end
end

%% --- Attach to table --------------------------------------------
T.Av_TensionPerDens       = Av_TensionPerDens;
T.Med_TensionPerDens      = Med_TensionPerDens;
T.Av_TensionPerDens_log10 = Av_TensionPerDens_log10;
T.N_TensionPerDens        = N_TensionPerDens;

end