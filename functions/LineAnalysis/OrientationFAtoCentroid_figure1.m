function OrientationFAtoCentroid_figure1(CELL, switchF, maskCC, centroidXY, centerFA, centerProximal, centerDistal)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   Still needs figure to be created


%% init
if ~switchF
    return
else
end

%% main

if exist('CELL', 'var')
    %%% create cell outline
    % get cell mask
    cellMask = CELL.mask;
    % dilate mask
    cellMaskDilate = imdilate(cellMask, strel('disk', 1));
    % remove inside of mask
    cellMaskOutline = cellMaskDilate - cellMask;
else
    cellMaskOutline = [];
end

for i = 1:maskCC.NumObjects

    if i == 1
        figure      
        title('Green stationary = centroid of all FAs | Green moving = centroid FA | blue = proximal | red = distal')
        hold on
    else        
    end

    % generate array to show FAs - this code below works nicely !
    verifyImage = zeros(maskCC.ImageSize);
    verifyVisualize = cast(verifyImage, 'uint8');
    verifyVisualize(maskCC.PixelIdxList{i}) = 255;
    verifyVisualize = cat(3, verifyVisualize, verifyVisualize, verifyVisualize);
    % show centroid
    centroidImage = zeros(maskCC.ImageSize);
    centroidImage(round(centroidXY(2)), round(centroidXY(1))) = 1;
    % centroidImage = imdilate(centroidImage, strel('disk', 5));
    tmp = verifyVisualize(:,:,3);
    tmp(centroidImage == 1) = 255;
    % show centroid of FA
    centroidImage(round(centerFA(i).Centroid(2)), round(centerFA(i).Centroid(1))) = 1;
    centroidImage = imdilate(centroidImage, strel('disk', 5));
    tmp(centroidImage == 1) = 255;
    verifyVisualize(:,:,2) = tmp;

    % show cell outline
    tmpr = verifyVisualize(:,:,1);
    tmpg = verifyVisualize(:,:,2);
    tmpb = verifyVisualize(:,:,3);
    tmpr(cellMaskOutline == 1) = 255;
    tmpg(cellMaskOutline == 1) = 255;
    tmpb(cellMaskOutline == 1) = 255;
    verifyVisualize(:,:,1) = tmpr;
    verifyVisualize(:,:,2) = tmpg;
    verifyVisualize(:,:,3) = tmpb;

    % show proximal end = BLUE
    centroidImage = zeros(maskCC.ImageSize);
    centroidImage(round(centerProximal(i, 1)), round(centerProximal(i, 2))) = 1;
    centroidImage = imdilate(centroidImage, strel('disk', 4));
    tmp = verifyVisualize(:,:,1);
    tmp(centroidImage == 1) = 255;
    verifyVisualize(:,:,3) = tmp;

    % show distal end = RED
    centroidImage = zeros(maskCC.ImageSize);
    centroidImage(round(centerDistal(i, 1)), round(centerDistal(i, 2))) = 1;
    centroidImage = imdilate(centroidImage, strel('disk', 4));
    tmp = verifyVisualize(:,:,1);
    tmp(centroidImage == 1) = 255;
    verifyVisualize(:,:,1) = tmp;

    % draw line
    % Bresenham's line algorithm is a line drawing algorithm
    [bres.x, bres.y] = bresenham(round(centerFA(i).Centroid(1)), round(centerFA(i).Centroid(2)), round(centroidXY(1)), round(centroidXY(2))); % input format: (x1,y1,x2,y2)
    tmpr = verifyVisualize(:,:,1);
    tmpg = verifyVisualize(:,:,2);
    tmpb = verifyVisualize(:,:,3);
    for k = 1:length(bres.x)
        % array format again → flip x and y
        tmpr(bres.y(k), bres.x(k)) = 101;  % **so simply flip x and y
        tmpg(bres.y(k), bres.x(k)) = 101;  % **so simply flip x and y
        tmpb(bres.y(k), bres.x(k)) = 101;  % **so simply flip x and y
    end
    verifyVisualize = cat(3, tmpr, tmpg, tmpb);


    pause(0.1)
    imshow(verifyVisualize)
    hold on
end
end

