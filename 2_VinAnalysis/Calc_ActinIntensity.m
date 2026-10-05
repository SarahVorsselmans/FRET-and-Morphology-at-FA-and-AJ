function T = Calc_ActinIntensity(T, options)
%CALC_ACTININTENSITY  Mask background in actin images and compute the average
%                     intensity of the remaining pixels, per image.
%
% USAGE
%   T_Act = Calc_ActinIntensity(T_Act);                          % fixed threshold 20
%   T_Act = Calc_ActinIntensity(T_Act, Threshold=22);            % fixed threshold 22
%   T_Act = Calc_ActinIntensity(T_Act, Method="isodata");        % automatic, per image
%   T_Act = Calc_ActinIntensity(T_Act, Method="otsu");
%
% OPTIONS
%   'Threshold'   : Fixed threshold (used when Method = "fixed"). Pixels with
%                   intensity >= Threshold are kept, pixels below become NaN.
%                   Default: 20.
%   'Method'      : "fixed" (default) | "isodata" | "otsu".
%                   isodata = iterative intermeans, close to Fiji's "Default"
%                   threshold (not identical: Fiji ignores the extreme bins).
%                   For the automatic methods, pixels > threshold are kept.
%   'ImageCol'    : Name of the image column. Default: 'ActinImage'.
%   'StoreMasked' : Logical. Store the masked image (background = NaN) in the
%                   column 'ActinMasked'. Default: true.
%
% OUTPUT COLUMNS ADDED TO T
%   AverageIntensity : mean of the pixels that remain after thresholding
%   ThresholdUsed    : threshold applied to that image
%   MaskedAreaFrac   : fraction of the image above threshold (0-1)
%   ActinMasked      : (optional) image with background set to NaN

arguments
    T table
    options.Threshold   (1,1) double = 20
    options.Method      string {mustBeMember(options.Method, ["fixed","isodata","otsu"])} = "fixed"
    options.ImageCol    string = "ActinImage"
    options.StoreMasked (1,1) logical = true
end

n          = height(T);
avgInt     = nan(n,1);
thrUsed    = nan(n,1);
areaFrac   = nan(n,1);
maskedImgs = cell(n,1);

for i = 1:n
    img = double(T.(options.ImageCol){i});

    switch options.Method
        case "fixed"
            thr  = options.Threshold;
            keep = img >= thr;
        case "isodata"
            thr  = isodataThreshold(img);
            keep = img > thr;
        case "otsu"
            thr  = otsuThresholdInt(img);
            keep = img > thr;
    end

    masked = img;
    masked(~keep) = NaN;

    avgInt(i)   = mean(masked, 'all', 'omitnan');   % NaN if nothing remains
    thrUsed(i)  = thr;
    areaFrac(i) = nnz(keep) / numel(img);
    maskedImgs{i} = masked;
end

T.AverageIntensity = avgInt;
T.ThresholdUsed    = thrUsed;
T.MaskedAreaFrac   = areaFrac;
if options.StoreMasked
    T.ActinMasked = maskedImgs;
end

fprintf('Actin intensity (%s threshold): mean threshold = %.1f (range %.1f-%.1f), mean area above threshold = %.1f%%\n', ...
    options.Method, mean(thrUsed), min(thrUsed), max(thrUsed), 100*mean(areaFrac));
end

% =========================================================================
% LOCAL HELPERS
% =========================================================================
function [counts, maxV] = intHistogram(img)
% Histogram with one bin per integer intensity level (0 ... max)
    vals = img(~isnan(img));
    maxV = round(max(vals));
    counts = histcounts(vals, -0.5:1:(maxV + 0.5));
end

function t = isodataThreshold(img)
% Iterative intermeans (IsoData) threshold on the integer histogram.
    [counts, maxV] = intHistogram(img);
    if isempty(counts) || maxV == 0, t = 0; return; end
    levels = 0:maxV;
    t = round(sum(levels .* counts) / sum(counts));   % start at the mean
    for it = 1:1000
        lo = levels <= t;
        hi = ~lo;
        if ~any(counts(lo)) || ~any(counts(hi)), break; end
        mLo = sum(levels(lo) .* counts(lo)) / sum(counts(lo));
        mHi = sum(levels(hi) .* counts(hi)) / sum(counts(hi));
        tNew = round((mLo + mHi) / 2);
        if tNew == t, break; end
        t = tNew;
    end
end

function t = otsuThresholdInt(img)
% Otsu threshold expressed in intensity units.
    [counts, maxV] = intHistogram(img);
    if isempty(counts) || maxV == 0, t = 0; return; end
    t = round(otsuthresh(counts) * (numel(counts) - 1));
end