function [output] = CellSegment(Int, threshFactor)
% CellSegment  - Segment whole-cell mask from fluorescence image
%   Uses morphological background estimation to suppress focal adhesions
%   and extract the cytosolic footprint.
%
% INPUTS:
%   Int          - raw intensity image (2D)
%   threshFactor - optional scaling factor (if NaN, default = 1.0)
%
% OUTPUT:
%   output       - binary mask of segmented cell

%% Handle threshFactor
if isnan(threshFactor)
    threshFactor = 1.0;  % neutral factor
end

%% Preprocessing
img = rescale(Int);         % scale intensities to [0,1]
img = medfilt2(img, [3 3]); % suppress salt & pepper noise

%% Step 1: Estimate cytosolic background by removing small bright structures (FAs)
SE = strel('disk', 20);     % radius defines "small" structures (tune for your data)
bg = imopen(img, SE);       % morphological opening = smooth background

%% Step 2: Threshold the background image
T = graythresh(bg);         % Otsu on smoothed background
bw = imbinarize(bg, T * threshFactor);

%% Step 3: Morphological cleanup
bw = imfill(bw, 'holes');             % fill gaps inside the cell
bw = bwareaopen(bw, 500);             % remove tiny specks
bw = imclose(bw, strel('disk', 10));  % smooth cell outline

%% Step 4: Keep largest connected region (the cell)
CC = bwconncomp(bw);
if CC.NumObjects > 0
    stats = regionprops(CC, 'Area');
    [~, idx] = max([stats.Area]);
    mask = false(size(bw));
    mask(CC.PixelIdxList{idx}) = true;
else
    mask = false(size(bw));
end

%% Output
output = mask;

end
