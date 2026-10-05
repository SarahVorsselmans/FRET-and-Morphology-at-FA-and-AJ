function [vectorFromCenterNorm] = OrientationDistance(show, vectorFromCenter, nBins, isBinNumber, binEdges)
%ORIENTATIONDISTANCE Summary of this function goes here
%   Detailed explanation goes here




%% Get euclidean distance of vectors from centroid to FA centroid
for i = 1:length(vectorFromCenter)
    vectorFromCenterNorm(i) = norm(vectorFromCenter{i});
end

%% find FA in each bin
for i = 1:ceil(nBins/2)
    k = find(isBinNumber == i);
    meanDistFromCenter(i) = mean(vectorFromCenterNorm(k));
end

%% Get angles
angleDifferenceHalf = rad2deg((abs(binEdges(1)) - abs(binEdges(2)) )/2) ;
for i = 1:floor(length(binEdges)/2)
    anglesDeg(i) = rad2deg(abs(binEdges(i))) - angleDifferenceHalf;
end

%% make figure
if show
    figure
    scatter(anglesDeg, meanDistFromCenter, 'filled')
    hold on
    plot(anglesDeg, meanDistFromCenter, '--')
    set ( gca, 'xdir', 'reverse' );
    xlim([0 180]);
    ylim([0 max(meanDistFromCenter)])
    ylabel(['Mean distance of FAs to center of cell' newline '(pixels)']);
    xlabel(['Absolute angle between FA and outward vector ' newline ' |' char(952) '| (deg)']);
else
end

%% parse
vectorFromCenterNorm = vectorFromCenterNorm';

end

