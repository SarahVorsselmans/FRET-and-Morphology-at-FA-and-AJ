function [] = intraFAanalysisPlot(lineFaIdxList, maskCC)
%INTRAFAANALYSISPLOT Summary of this function goes here

%     Copyright (C) 2024 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2024 Laurens Kimps
%

%% parse

% mask
maskIdxList{:,1} = cat(1, maskCC.PixelIdxList{1:end});

% create array from cell
maskIdxList = cell2mat(maskIdxList);

%% plot function • moving static

% create static image
imgStatic = cast(false(maskCC.ImageSize), 'uint8');
imgStaticRed = imgStatic;
imgStaticBlue = imgStatic;
imgStatic = repmat(imgStatic, [1 1 3]);
imgStaticRed(lineFaIdxList) = 1; % lineFA = red
imgStaticBlue(maskIdxList) = 1; % mask = blue
imgStaticBlue(imgStaticRed == 1) = 0;
imgStatic(:,:,3) = imgStaticRed;
imgStatic(:,:,1) = imgStaticBlue;
imgStatic = rescale(imgStatic);

% create figure
f = figure;
f.Position(1:2) = [150 150];
f.Position(3:4) = [750 800];
f.Name = "FA line analysis • static FAs";
imshow(imgStatic)
title('There should be no red pixels in the image')


%% plot function • moving segments

% create image
imgMove = cast(false(maskCC.ImageSize), 'uint8');
imgMove(lineFaIdxList) = 100;

% create figure
f = figure;
f.Position(1:2) = [100 100];
f.Position(3:4) = [750 800];
f.Name = "FA line analysis • moving segments";
imshow(imgMove)

% titleText = join(["•   FA # 1/" num2str(maskCC.NumObjects) ...
%     "  •   segment 1/" num2str(length(lineFAnonZero(i_FA).segmentPixels{1}))]); 
% title(titleText)

% loop through FA segments
for i_FA = 1:maskCC.NumObjects
    for i_line = 1:length(lineFAnonZero(i_FA).segmentPixels)
        imgMove(lineFAnonZero(i_FA).segmentPixels{i_line}) = 255;
        imshow(imgMove)
%         titleText = join(["•   FA # " num2str(i_FA) "/" num2str(maskCC.NumObjects) ...
%             "  •   segment # " num2str(i_line) "/" num2str(length(lineFAnonZero(i_FA).segmentPixels{i_line}))]); 
%         title(titleText)
        pause(0.010)
    end
end



end