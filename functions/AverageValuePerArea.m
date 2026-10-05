function [averageIndexPerArea, areaInPixels, connComp] = AverageValuePerArea(data, mask, delIdx, flagShowRemoved, cnt)
%FCN_LARGESTAREA produces the average FRET index per area, the size of each
%area in pixels and gives the connected component as a structure

%   <data> is the fret index array
%   <mask> is the binary array of the masked fret index

% set NaNs to zero
data(isnan(data)) = 0;

% get separate areas
connComp = bwconncomp(mask, 8);
areaInPixels = regionprops(connComp, 'Area');
statsMeanIntensity = regionprops(connComp, data,'MeanIntensity');
[averageIndexPerArea] = deal([statsMeanIntensity.MeanIntensity]);
averageIndexPerArea = averageIndexPerArea';

oldRegions = zeros([connComp.ImageSize]);
imgLbl = zeros([connComp.ImageSize]);
for i = 1:connComp.NumObjects
    oldRegions(connComp.PixelIdxList{i}) = true;
    imgLbl(connComp.PixelIdxList{i}) = i;
end


connComp.PixelIdxList(delIdx) = [];
connComp.NumObjects = length(connComp.PixelIdxList) ;
areaInPixels(delIdx) = [];
averageIndexPerArea(delIdx) = [];

if flagShowRemoved & cnt & any(delIdx)
    newRegions = zeros([connComp.ImageSize]);
    for i = 1:connComp.NumObjects
        newRegions(connComp.PixelIdxList{i}) = true;
    end

    imFused = imfuse(oldRegions, newRegions);
    figure("Name", 'Deleted FAs')
    imshow(imFused)
    title('You removed these FAs (green)')
    cnt = 1;
else
end

% make figure that tells you which FA is which with the data cursor
if cnt == 0
    figure("Name", 'FA index finder')
    imshow(imgLbl)
    title('If you hover over the FA you can see its index = FA label (use "data tips" tool)')

    cnt = 1;
else
end

end
