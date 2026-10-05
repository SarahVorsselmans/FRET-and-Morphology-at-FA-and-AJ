function [lineFAnonZero, lineFaIdxList] = intraFAanalysisBase(maskCC)
%INTRAFAANALYSISBASE Analyze within binary regions

%     Copyright (C) 2024 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2024 Laurens Kimps
% 
%     This code makes use of the bresenham.m function
%     Copyright (c) 2010, Aaron  Wetzler
%     All rights reserved. (see license.txt)
%     
%     <maskCC>: bwconncomp of binary mask
%     <lineFA>: a structure with a list of pixel indices for each region
%     for each FA

%% Main

disp("• Partitioning FAs ...")

% find regionprops
stats = regionprops(maskCC, "BoundingBox", "MajorAxisLength", "Orientation", "Centroid");

% loop over number of binary elements
for i = 1:maskCC.NumObjects
    
    disp([" - Partitioning FA " num2str(i) "/" num2str(maskCC.NumObjects) " into segments ..."])

    % for each element (FA):
    
    % reduce matrix to the data within bouding box
    region1 = false(maskCC.ImageSize);
    region1(maskCC.PixelIdxList{i}) = 1;
    
    % get bouding box coordinates (for length perpindicular line)
    bres_x = floor(stats(i).BoundingBox(1)); % matlab image convention: x →
    horz = ceil(stats(i).BoundingBox(3)); % height -> flip
    bres_y = floor(stats(i).BoundingBox(2)); % matlab image convention: y ↓
    vert = ceil(stats(i).BoundingBox(4)); % width -> flip
    
    % image format
    x_var = (bres_x+1:bres_x+horz)+0;
    y_var = (bres_y+1:bres_y+vert)+0;

    % Matlab array →    Matlab image
    % x: vert down →    x: hor right
    % y: hor right →    y: vert down
    % **so simply flip x and y 
    
    % Matlab array format: x ↓ and y →
    region2 = region1(y_var, x_var); % array: x = vert down | y = hor right
%     dipshow(region2)
%     figure, imshow(region1)

    % define major axis (what if region is a circle?)
    
    % matlab image convention: x → and y ↓
    xj = stats(i).Centroid(1) + [-1 1] * (stats(i).MajorAxisLength + 5) * cosd(stats(i).Orientation) / 2;
    yj = stats(i).Centroid(2) - [-1 1] * (stats(i).MajorAxisLength + 5) * sind(stats(i).Orientation) / 2;
    % → start from centroid
    % → create vector going in plus and minus direction
    % → multiply by length divided by 2 and corrected for the angle
    % Added 5 pixel(s) to the MajorAxisLength to avoid undersampling

    % include folders to path
    addpath("bin\functions\community\bresenham\");
    
    % Bresenham's line algorithm is a line drawing algorithm
    [bres.x, bres.y] = bresenham(xj(1), yj(1), xj(2), yj(2)); % input format: (x1,y1,x2,y2)
    
    % make logical → uint8 (for labeling)
    lineimg = cast(region1, 'uint8'); 
    % label line as value 2
    for k = 1:length(bres.x)
        % array format again → flip x and y
        lineimg(bres.y(k), bres.x(k)) = 2;  % **so simply flip x and y
    end

    % crop bounding box
    %lineimg = lineimg(y_var, x_var); % **so simply flip x and y

    %     dipshow(lineimg);
            
    % define longest possible line in bounding box
    diag = sqrt(horz^2 + vert^2);
    
    % figure for testing
%     figure
%     hold on    

    % Define lines image
    perpLines = false(maskCC.ImageSize);
    sampleRegion = false(maskCC.ImageSize);
    previousSampleRegion = false(maskCC.ImageSize);;
    samplShow = cast(region1, 'uint8');
    perpBackup = perpLines;

    % walk along major axis & create perpindicular axis

    for m = 1:length(bres.x)

        % → 'centroid' is any point along line
        % →  multiply for vector in both directions
        % change the angle 90°
        xj = bres.x(m)  +  [-1 1] * diag*cosd(stats(i).Orientation + 90 )/2;
        yj = bres.y(m)  -  [-1 1] * diag*sind(stats(i).Orientation + 90 )/2;

        % Bresenham's line algorithm is a line drawing algorithm
        [perp.x, perp.y] = bresenham(xj(1), yj(1), xj(2), yj(2)); % input format: (x1,y1,x2,y2)

        % label perp line as value 3
        for k = 1:length(perp.x)
            % array format again → flip x and y
            lineimg(perp.y(k), perp.x(k)) = 3;  % **so simply flip x and y
            perpLines(perp.y(k), perp.x(k)) = true;
        end
        
%         figure for testing
%         imshow(label2rgb(lineimg,'jet','k','shuffle'))
%         pause(0.1)

        % clean perpindicular sampling region
        perpLines(perpBackup == 1) = 1; % needed to fill all holes
        perpLines = imfill(perpLines, 'holes');
        perpBackup = perpLines;
        perpLines(region1 == 0) = 0;        
        sampleRegion(perpLines == 1) = 1;
        ExtraAddedPixels = sampleRegion - previousSampleRegion;
        previousSampleRegion = sampleRegion;

        % generate indexable FA line pixel list-array in structure
        lineFA(i).segmentPixels{m} = find(ExtraAddedPixels == 1);
        
%         % figure for testing
%         samplShow(perpLines == 1) = 2;
%         imshow(label2rgb(samplShow,'jet','k','shuffle'))
%         pause(0.1)
        
%         figure for testing
%         imshow(sampleRegion);
%         pause(0.1)

    end

end

disp("• Partitioning FAs done")

%% generate mask from indices (+ remove empty cell entries)
% to verify that all pixels of maskCC have been alocated  to a line segment

% initialize variables
tmp_01 = {};
lineFAnonZero = lineFA;
maskIdxList = {};

% loop through FAs
for i = 1:maskCC.NumObjects

    % for lineFA
    tmp_03 = lineFA(i).segmentPixels; % create new variable for readability
    tmp_03 = tmp_03(~cellfun('isempty', tmp_03)); % remove [] entries
    lineFAnonZero(i).segmentPixels = tmp_03; % update lineFA
    tmp_01{i,1} = cat(1, lineFAnonZero(i).segmentPixels{1:end});
    % cat also removes empty entries from the list

end

% for mask
maskIdxList{:,1} = cat(1, maskCC.PixelIdxList{1:end});

% create array from cell
lineFaIdxList = cell2mat(tmp_01);
maskIdxList = cell2mat(maskIdxList);


end

%         % Alternative method might not work better 
%         bwPerpLine = false(maskCC.ImageSize);
% 
%         for k = 1:length(perp.x)
%             % array format again → flip x and y
%             bwPerpLine(perp.y(k), perp.x(k)) = true;  % **so simply flip x and y
%         end
% 
%         for k = 1:length(perp.x)
%             SE = strel("disk", 3);
%             bwPerpLine = imdilate(bwPerpLine, SE);
% 
%             %             % figure for testing
%             %             imshow(bwPerpLine)
%             %             pause(0.1)
% 
%             % figure for testing
%             lineimg(bwPerpLine == 1) = 5;
%             lineimg(region1 == 1) = 1;
%             imshow(label2rgb(lineimg,'jet','k','shuffle'))
%             pause(0.1)
%             % End alternative method
% 
%         end