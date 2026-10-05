function OrientationFAtoCentroid_figure2(testMode, binEdges, maskCC, isBinNumber)
%ORIENTATIONFATOCENTROID_FIGURE2 Summary of this function goes here
%   Detailed explanation goes here

%% run?

if ~testMode
    return
else
end

%% start

figure

% loop through bins
for b3 = 1:floor(length(binEdges)/2)
    k = find(isBinNumber == b3);

    limitLow1 = binEdges(b3);
    limitHigh2 = binEdges(end + 1 - b3);

    limitHigh1 = binEdges(b3 + 1);
    limitLow2 = binEdges(end - b3);

    limitString = ['Angle intervals: [' num2str(limitLow1) ' ' num2str(limitHigh1) '] or [' num2str(limitLow2) ' ' num2str(limitHigh2) ']' ];

    if b3 == 1
        % generate array to show FAs
        verifyImage = zeros(maskCC.ImageSize);
        verifyVisualize = cast(verifyImage, 'uint8');
        for b5 = 1:maskCC.NumObjects
            verifyVisualize(maskCC.PixelIdxList{b5}) = 50;
        end
        verifyVisualize = cat(3, verifyVisualize, verifyVisualize, verifyVisualize);
    else
    end

    tmp = verifyVisualize(:,:,1);

    for b4 = 1:length(k)
        tmp_idx = maskCC.PixelIdxList{k(b4)};
        tmp(tmp_idx) = 255;
    end
    verifyVisualize = cat(3, tmp, verifyVisualize(:,:,2), verifyVisualize(:,:,2));

    imshow(verifyVisualize)
    title(limitString)

    pause(0.5)

    verifyVisualize(verifyVisualize == 255) = 80;

end
end

