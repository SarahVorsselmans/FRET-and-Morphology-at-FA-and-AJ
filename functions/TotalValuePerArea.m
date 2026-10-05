function [totalIntensity] = TotalValuePerArea(data, cc)
%FCN_LARGESTAREA produces the average FRET index per area, the size of each
%area in pixels and gives the connected component as a structure

%   <data> is the fret index array
%   <mask> is the binary array of the masked fret index

% set NaNs to zero
data(isnan(data)) = 0;

% get total value per area
for i = 1:cc.NumObjects
    totalIntensity(i) = sum(data(cc.PixelIdxList{i}));
end
