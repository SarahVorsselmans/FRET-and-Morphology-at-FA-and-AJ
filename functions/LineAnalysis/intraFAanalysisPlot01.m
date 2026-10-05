function intraFAanalysisPlot01(switchF, maskCC, lineIdx, centroid)
%INTRAFAANALYSISPLOT01 Summary of this function goes here

%     Copyright (C) 2024 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2024 Laurens Kimps

%% run ?
if ~switchF
    return
else
end

%% start

% visualize centroid in array
centroidImage = zeros(maskCC.ImageSize);
centroidImage(round(centroid(2)), round(centroid(1))) = 1;
centroidImage = imdilate(centroidImage, strel('disk', 5));

% prepare variables
verifyImage = false(maskCC.ImageSize);
verifyImageLine = verifyImage;
verifyImageMask = verifyImage;

C2 = [lineIdx{:}];
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
disp(['     - Overlap between FA mask and FA partitions = ' num2str(U) ' (should be 1)'])

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

for i = [3 2 1]
    title([num2str(i)])
    pause(1)
end

title('Blue dot = centroid of FAs')

% for w = 1:length({lineFA(1:end).segmentPixels})
for w1 = 1:max(cellfun(@length, lineIdx)) % w1 = line segment of FA
    % for longest FA line section

    for w2 = 1:maskCC.NumObjects % w2 = loop FAs

        if w1 > length(lineIdx{w2}) % skip current FA
            continue
        else
            tmp_vis_1 = lineIdx{w2}(w1); % select FA
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
    title('Blue dot = centroid of FAs')
    hold on

end

end

