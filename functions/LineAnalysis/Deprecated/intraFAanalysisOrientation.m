function [lineFAnonZero] = intraFAanalysisOrientation(lineFAnonZero, cellMask)
%INTRAFAANALYSISORIENTATION Summary of this function goes here
%   Detailed explanation goes here

%     Copyright (C) 2024 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2024 Laurens Kimps
%

%% Main

% Get centroid location
% centroid = ...

% image size
ImageSize = maskCC.ImageSize;

disp("• Checking FA orientation with respect to centroid location")

% loop through FAs
for i_FA = 1:length(lineFAnonZero.segmentPixels)
    
    % find first and last pixel index of each FA
    firstPixel(i_FA) = lineFAnonZero(i_FA).segmentPixels{1};
    lastPixel(i_FA) = lineFAnonZero(i_FA).segmentPixels{end};

    xy_first = ind2sub(ImageSize, firstPixel(i_FA));
    xy_last = ind2sub(ImageSize, lastPixel(i_FA));

    distFirst = sqrt((xy_first(1)-centroid(1)).^2+(xy_first(2)-centroid(2)).^2);
    distLast = sqrt((xy_last(1)-centroid(1)).^2+(xy_last(2)-centroid(2)).^2);

    % arrange from distal to proximal
    if distFirst >= distLast
        % keep distal point first in list: do nothing
        flipTracker(i_FA) = 1;
        disp([" - FA #" num2str(i_FA) "/" num2str(length(lineFAnonZero.segmentPixels)) " → ok"])

    elseif distFirst <= distLast
        % flip to have distal point first in list: flip
        lineFAnonZero(i_FA).segmentPixels = flip(lineFAnonZero(i_FA).segmentPixels);
        flipTracker(i_FA) = -1;
        disp([" - FA #" num2str(i_FA) "/" num2str(length(lineFAnonZero.segmentPixels)) " → corrected orientation"])

    else
        flipTracker(i_FA) = 0;
        disp([" - FA #" num2str(i_FA) "/" num2str(length(lineFAnonZero.segmentPixels)) " → Tangential to centroid"])

    end

%     <flipTracker> shows the initial orientation of the FA pixel Idx list

disp("• FA orientation check: done")
end

