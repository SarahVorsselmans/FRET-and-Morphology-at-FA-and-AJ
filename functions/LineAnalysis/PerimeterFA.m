function [PerimeterIndex] = PerimeterFA(CELL, fa, outputPath)
%PERIMETERFA Summary of this function goes here
%   Detailed explanation goes here


%% Overall perimeter

% get mask
mask = CELL.mask;
% get conncomp
CC1 = bwconncomp(mask);
% get Minimum Feret Diameter
stats = regionprops(CC1, 'MinFeretProperties');
% define factor for perimeter width
perimFactor = 0.09;
% calculate perimeter width
perimWidth = perimFactor * stats.MinFeretDiameter;
% get smaller inside area
tmpMaskPerim = imerode(mask, strel("disk", round(perimWidth)));
% substract smaller area from bigger area = perimeter
maskPerim = cast((mask - tmpMaskPerim), 'logical');

peri = bwconncomp(maskPerim);

if 0
    dipshow(maskPerim)
else
end


%% Main

% init vars
binaryCheck = [];
belongsToRegion = [];

faimg = zeros(fa.ImageSize);
faperi = zeros(fa.ImageSize);

% loop over FA regions
for i = 1:length(fa.PixelIdxList)

    testa = zeros(fa.ImageSize);
    testb = zeros(fa.ImageSize);


    flimIdxList = fa.PixelIdxList{i};
    testa(flimIdxList) = 1;
    faimg(flimIdxList) = 1;

    regionContainer = {};
    binaryCheckRegionNum = [];

    % loop over field regions
    for ii = 1:peri.NumObjects
        tfmIdxList = peri.PixelIdxList{ii};
        testb(tfmIdxList) = 1;
        regionContainer{ii} = [ismember(flimIdxList, tfmIdxList)];
        binaryCheckRegionNum(ii) = any(regionContainer{ii} == 1);
    end

    % check for ones (and set to zero if none found in tmp)
    if any(binaryCheckRegionNum(:) == 1) == 0
        belongsToRegion = [belongsToRegion 0];
    else
        % store disp / tract field label in FA conncomp list
        [~,idx] = ismember(1, binaryCheckRegionNum(:) );
        belongsToRegion = [belongsToRegion idx];
        faperi(flimIdxList) = 255;
    end
    binaryCheck = [binaryCheck any([binaryCheckRegionNum(:)] == 1)];
end
PerimeterIndex = belongsToRegion';

%% create image

empty = zeros(fa.ImageSize);

bm = empty;
bm(maskPerim == 1) = 255;

r = empty;
g = empty;
b = empty;

b(faimg == 1) = 255;
g(faimg == 1) = 100;

r(bm == 255) = 255; 
g(bm == 255) = 50; 
b(bm == 255) = 50;  


r(faperi == 255) = 255; 
g(faperi == 255) = 255; 
b(faperi == 255) = 255;  

img = cat(3, r, g, b);

img = cast(img, 'uint8');
figure, imshow(img)
title('red: perimeter, white: FAs on perimeter, blue: cytoplasmic FAs')

saveas(gca, [outputPath filesep 'Perimeter.png'])

end

