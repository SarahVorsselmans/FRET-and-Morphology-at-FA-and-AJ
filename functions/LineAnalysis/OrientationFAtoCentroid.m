function [vectorFromCenter, nBins, isBinNumber, binEdges, FAcoordinates, angleDeg] = OrientationFAtoCentroid(maskCC, idxList, centroidXY)
%ORIENTATIONFATOCENTROID Summary of this function goes here
%     Copyright (C) 2024 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2024 Laurens Kimps
%
%
%     <maskCC>: bwconncomp of binary mask
%     <lineFA>: a structure with a list of pixel indices for each region
%     for each FA

%%
testMode = 0;

%% main

% calc vector from centroid to FA center
% calc vector from FA proximal to distal end = bwconncomp stats no!
% calc difference between these vectors
% make polar plot


centerFA = regionprops(maskCC, 'centroid');
ImageSize = maskCC.ImageSize;

% loop through FAs
for i = 1:length(idxList)

    % Reminder:
    % Matlab array →    Matlab image
    % x: vert down →    x: hor right
    % y: hor right →    y: vert down
    % **so simply flip x and y

    % calc vector from centroid to FA center
    vectorFromCenter{i} = flip(centerFA(i).Centroid) - flip(centroidXY);

    % calc vector from FA proximal to distal end = bwconncomp stats no!
    centerProximal_tmp = idxList{i}{1}; % get center pixel
    centerDistal_tmp = idxList{i}{end}; % get center pixel

    [centerProximalX, centerProximalY] = (ind2sub(ImageSize, centerProximal_tmp)); % proximal X Y
    [centerDistalX, centerDistalY] = (ind2sub(ImageSize, centerDistal_tmp)); % distal X Y

    centerProximal(i, 1:2) = round([mean(centerProximalX) mean(centerProximalY)]); % get centroid
    centerDistal(i, 1:2) = round([mean(centerDistalX) mean(centerDistalY)]); % get centroid

    vectorProx2Dist{i} = centerDistal(i, 1:2) - centerProximal(i, 1:2); % get vector pointing from proximal to distal end of FA

    %     % parse
    %     x1 = vectorFromCenter{i}(1);
    %     x2 = vectorFromCenter{i}(2);
    %     y1 = vectorProx2Dist{i}(1);
    %     y2 = vectorProx2Dist{i}(2);

    %     % find angle between vector from centroid and from proximal to distal end
    %     angleDeg(i) = atan2d( (  (vectorFromCenter{i}(1) * vectorProx2Dist{i}(2)) - (vectorFromCenter{i}(2) * vectorProx2Dist{i}(1))  ) / (dot(vectorFromCenter{i}, vectorProx2Dist{i})) );
    %     angleRad(i) = atan2( (  (vectorFromCenter{i}(1) * vectorProx2Dist{i}(2)) - (vectorFromCenter{i}(2) * vectorProx2Dist{i}(1))  ) / (dot(vectorFromCenter{i}, vectorProx2Dist{i})) );

    % angleDeg(i) = atan2d(vectorFromCenter{i}, vectorProx2Dist{i});
    % angleRad(i) = atan2(vectorFromCenter{i}, vectorProx2Dist{i});
    %

    angleDeg(i) = atan2d( (  (vectorFromCenter{i}(1) * vectorProx2Dist{i}(2)) - (vectorFromCenter{i}(2) * vectorProx2Dist{i}(1))  ), (dot(vectorFromCenter{i}, vectorProx2Dist{i})) );
    angleRad(i) = atan2( (  (vectorFromCenter{i}(1) * vectorProx2Dist{i}(2)) - (vectorFromCenter{i}(2) * vectorProx2Dist{i}(1))  ), (dot(vectorFromCenter{i}, vectorProx2Dist{i})) );



    % %     % find angle between vector from centroid and HORIZONTAL Y VECTOR
    % %     vectorHorz = [0 1];
    %     angleDegHorz(i) = acosd( (dot(vectorFromCenter{i}, vectorHorz)) / (norm(vectorFromCenter{i})*norm(vectorHorz)) );
    %     angleRadHorz(i) = acos( (dot(vectorFromCenter{i}, vectorHorz)) / (norm(vectorFromCenter{i})*norm(vectorHorz)) );
    %
    %     angleDeg(i) = angleDeg(i) + angleDegHorz(i);
    %     angleRad(i) = angleRad(i) + angleRadHorz(i);
    %

end

    %% Validation | visualize centers
    if testMode      
        OrientationFAtoCentroid_figure1(1, maskCC, centroidXY, centerFA, centerProximal, centerDistal)
    else
    end

%% generate polarhistogram

nBins = 4*12+1; % 12: this generates 15 bins between +90° and -90° of 15° each

figure
h = polarhistogram(angleRad, nBins,'BinLimits', [-pi pi]);
title('Angle between FA and vector pointing from overall FA centroid to FA')
% polarhistogram(angleRad, 'NumBins', 30)


% validate by plotting the FAs based on their angle in a specific color



if mod(length(h.BinEdges), 2) == 0
    edges = [0 (h.BinEdges((length(h.BinEdges) /2) +1 :end)) ];
elseif mod(length(h.BinEdges), 2) ~= 0
    edges = h.BinEdges(ceil(length(h.BinEdges) /2 ) : end);
end


% Find which FA corresponds to which bin
for b1 = 1:length(idxList)
    for b2 = 1:length(edges - 1)
        if abs(angleRad(b1)) > edges(b2) && abs(angleRad(b1)) < edges(b2+1)
            isBinNumber(b1) = b2;
            isBinLimits{b1} = [edges(b2), edges(b2+1) ];
        else
        end
    end
end





%% Validation | visualize angles

if testMode    
    OrientationFAtoCentroid_figure2(testMode, h.BinEdges, maskCC, isBinNumber);
else
end

%% outputs
binEdges = h.BinEdges;
FAcoordinates.FaCentroids = centerFA;
FAcoordinates.ProximalVec = centerProximal;
FAcoordinates.DistalVec = centerDistal;
FAcoordinates.vectorProximal2Distal = vectorProx2Dist;
end








