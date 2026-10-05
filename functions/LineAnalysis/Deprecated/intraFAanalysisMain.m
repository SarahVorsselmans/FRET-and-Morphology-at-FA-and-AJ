function [lineFaIdxList, centroidXY] = intraFAanalysisMain(maskCC, user, CELL)
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

%% Testing

testMode = 1;

%% Init

% include folders to path
addpath("functions\community\bresenham\");

%% Main: partition FAs

disp('• Partitioning FAs ...')

% find regionprops
stats = regionprops(maskCC, "BoundingBox", "MajorAxisLength", "Orientation", "Centroid");

% loop over number of binary elements
for i = 1:maskCC.NumObjects
    % i = a focal adhesion

    disp(['     - Partitioning FA ' num2str(i) '/' num2str(maskCC.NumObjects) ' into line segments ...'])

    % for each element (FA ~ i):

    % reduce data in matrix to single FA (or region)
    region1 = false(maskCC.ImageSize);
    region1(maskCC.PixelIdxList{i}) = 1;

    % get bouding box coordinates (for length perpindicular line)
    bres_x = floor(stats(i).BoundingBox(1)); % matlab image convention: x →
    horz = ceil(stats(i).BoundingBox(3)); % height -> flip
    bres_y = floor(stats(i).BoundingBox(2)); % matlab image convention: y ↓
    vert = ceil(stats(i).BoundingBox(4)); % width -> flip

    % image format
    x_var = (bres_x+1:bres_x+horz);
    y_var = (bres_y+1:bres_y+vert);

    % Reminder:
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
    % xj line coordinate 1
    xj = stats(i).Centroid(1) + [-1 1] * (stats(i).MajorAxisLength + 5) * cosd(stats(i).Orientation) / 2;
    % yj line coordinate 2
    yj = stats(i).Centroid(2) - [-1 1] * (stats(i).MajorAxisLength + 5) * sind(stats(i).Orientation) / 2;
    % → start from centroid
    % → create vector going in plus and minus direction
    % → multiply by length divided by 2 and corrected for the angle
    % Added 5 pixel(s) to the MajorAxisLength to avoid undersampling

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

    %     if testMode
    %         % figure for testing
    %         figure
    %         hold on
    %     else
    %     end

    % Define lines image
    perpLines = false(maskCC.ImageSize);
    sampleRegion = false(maskCC.ImageSize);
    previousSampleRegion = false(maskCC.ImageSize);;
    samplShow = cast(region1, 'uint8');
    perpBackup = perpLines;

    % walk along major axis of FA & create perpindicular line
    for m = 1:length(bres.x)
        % m = centroid = point on the line
        % → 'centroid' is now any point along line
        % →  multiply for vector in both directions
        % change the angle 90°
        xj = bres.x(m)  +  [-1 1] * diag*cosd(stats(i).Orientation + 90 )/2;
        yj = bres.y(m)  -  [-1 1] * diag*sind(stats(i).Orientation + 90 )/2;

        % Bresenham's line algorithm is a line drawing algorithm
        [perp.x, perp.y] = bresenham(xj(1), yj(1), xj(2), yj(2)); % input format: (x1,y1,x2,y2)

        % label perp line as value 3
        for k = 1:length(perp.x)
            % k = index to point along the perpindicular line
            % array format again → flip x and y
            lineimg(perp.y(k), perp.x(k)) = 3;  % **so simply flip x and y
            perpLines(perp.y(k), perp.x(k)) = true;
        end

        %         if testMode
        %             figure for testing
        %             imshow(label2rgb(lineimg,'jet','k','shuffle'))
        %             pause(0.1)
        %         else
        %         end

        % clean perpindicular sampling region
        perpLines(perpBackup == 1) = 1; % needed to fill all holes
        perpLines = imfill(perpLines, 'holes');
        perpBackup = perpLines;
        perpLines(region1 == 0) = 0; % exclude datapoints outside of FA (or region)
        sampleRegion(perpLines == 1) = 1;
        ExtraAddedPixels = sampleRegion - previousSampleRegion;
        previousSampleRegion = sampleRegion;

        % generate list-array in struct with indexes of perp-FA-line pixels
        lineFA(i).segmentPixels{m} = find(ExtraAddedPixels == 1);

        %         if testMode
        %                     figure for testing
        %                     samplShow(perpLines == 1) = 2;
        %                     imshow(label2rgb(samplShow,'jet','k','shuffle'))
        %                     pause(0.1)
        %
        %                     figure for testing
        %                     imshow(sampleRegion);
        %                     pause(0.1)
        %         else
        %         end

    end
end

disp('• Partitioning FAs done')

%% Prepare FA partitions

% go from struct to array
C1 = {lineFA(1:end).segmentPixels};

% remove nonzero cell entries
for i = 1:length(C1)
    C1nz{i} = C1{i};
    C1nz{i} = C1nz{i}(~cellfun('isempty', C1nz{i})); % remove [] entries
end

%% Order FA partitions from proximal to distal
% get centroids
centroidsAll = {stats(1:end).Centroid};
centroidsAll = centroidsAll(:);
centroidsAll = cat(1, centroidsAll{1:end});
centroidX = mean(centroidsAll(1:end, 1));
centroidY = mean(centroidsAll(1:end, 2));

centroidXY = [centroidX centroidY];

% Reminder:
% Matlab array →    Matlab image
% x: vert down →    x: hor right
% y: hor right →    y: vert down
% **so simply flip x and y

% visualize centroid in array
centroidImage = zeros(maskCC.ImageSize);
centroidImage(round(centroidY), round(centroidX)) = 1;
centroidImage = imdilate(centroidImage, strel('disk', 5));

% change orientation if needed
[C1nz] = intraFAanalysisOrientation(C1nz, [centroidY centroidX], maskCC);

%% validation | generate mask from indices
% to verify that all pixels of maskCC have been alocated  to a line segment

if testMode

    % prepare variables
    verifyImage = false(maskCC.ImageSize);
    verifyImageLine = verifyImage;
    verifyImageMask = verifyImage;

    C2 = [C1nz{:}];
    % C2 = C2(~cellfun('isempty', C2)); % remove [] entries
    C2 = cat(1, C2{1, 1:end});

    verifyImageLine(C2) = 1;

    % go from struct to array
    C3 = {maskCC(1:end).PixelIdxList};
    C3 = [C3{:}];
    C3 = cat(1, C3{1, 1:end});

    verifyImageMask(C3) = 1;

    % figure
    % joinchannels('rgb', verifyImageLine, verifyImageMask)
    U = isequal(verifyImageLine, verifyImageMask);
    % U should be '1', nothing else
    disp('• Verifying FA partitioning:')
    disp(['  - Overlap between FA mask and FA partitions = ' num2str(U) ' (should be 1)'])

    %% Visualize the line segments on the FAs

    figure
    hold on

    % generate array to show FAs
    verifyVisualize = cast(verifyImage, 'uint8');
    verifyVisualize(verifyImageMask == 1) = 50;
    verifyVisualize = cat(3, verifyVisualize, verifyVisualize, verifyVisualize);
    % show centroid
    tmp = verifyVisualize(:,:,3);
    tmp(centroidImage == 1) = 255;
    verifyVisualize(:,:,3) = tmp;
    imshow(verifyVisualize)

    % for w = 1:length({lineFA(1:end).segmentPixels})
    for w1 = 1:max(cellfun(@length, C1nz)) % w1 = line segment of FA
        % for longest FA line section

        for w2 = 1:length({lineFA(1:end).segmentPixels}) % w2 = loop FAs

            if w1 > length(C1nz{w2}) % skip current FA
                continue
            else
                tmp_vis_1 = C1nz{w2}(w1); % select FA
                tmp_vis_1 = cat(1, tmp_vis_1{1, 1:end});
                %             verifyVisualize(tmp_vis_1) = 255;
                r =  verifyVisualize(:,:,1);
                g =  verifyVisualize(:,:,2);
                b =  verifyVisualize(:,:,3);

                r(tmp_vis_1) = 255;
                g(tmp_vis_1) = 255;
                b(tmp_vis_1) = 255;

                verifyVisualize = cat(3, r, g, b);
            end
        end
        if w1 <= 5
            pause(1)
        else
        end
        pause(0.1)
        imshow(verifyVisualize)
        hold on

    end

else
end

%% parse output

lineFaIdxList = C1nz;


