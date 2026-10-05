function [centroid, mask, maskAreaPixels] = CentroidAreaOfCellFinder(image)

%     Copyright (C) 2024 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2024 Laurens Kimps


%     This function assumes the TFM intensity image as input <image> (more
%     noisy and a projection of multiple stacks → less gaps between FAs)
%     I chose a threshold value of intensity <= 10 as the background signal

%% Init


%% main - segment image

% apply medfiltering 
J = medfilt2(image,[15 15]);
% duplicate variable
J2 = J;
% remove background
J2(J2 < 5) = 0;
% Otsu threshold
bw = imbinarize(J2);
% generate conncomp
CC = bwconncomp(bw); %make a list of all the objects in your mask
% find areas of objects
stats = regionprops(CC,'area');%gives you the area of each objective
% create 1D array from cell array
[areas] = deal([stats.Area]); 
% find index of largest object 
[~, indx] = max(areas); %just tells you which object has the maximum area, and you get the index information of this object
% create empty image
bw2 = false(size(bw));%you create an empty mask (everything is false or Zero)
% create binary mask of largest object
bw2(CC.PixelIdxList{indx}) = true;%then you make only the object with maximum area 1 ('true')
% create strel disk
SE = strel('disk',0);%creates a shape, but now its set to Zero 
% dilate image
bw2 = imdilate(bw2,SE);%applies the shape everywhere where the index is 1
% fill holes
bw2 = imfill(bw2,'holes');%any closed shape with hole will be filled
% erode again (to end up with same size)
mask = imerode(bw2,SE);             % <----------------------- output <mask>


if 0
    dipshow(mask)
else
end
%% main - centroid/area

% get conncomp of mask
CCMask = bwconncomp(bw2);
% get properties
statsMask = regionprops(CCMask,'Area', 'Centroid');
% find centroid
statsCentroid = statsMask.Centroid;
% round centroid values (actual pixel coordinates)
centroid = round(statsCentroid);    % <----------------------- output <centroid>
% find Area (pixels)
maskAreaPixels = statsMask.Area;    % <----------------------- output <maskAreaPixels>


%% End



end

