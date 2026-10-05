function [angleDeg, orientationAvgTFM] = OrientationFAtractions(showFigures, maskCC, listIdx, targetX, targetY, FA)
%ORIENTATIONFATRACTIONS Summary of this function goes here
%   Detailed explanation goes here

%% Parse

vectorProx2Distal = FA.Coordinates.vectorProximal2Distal;

tmp = struct2cell(FA.Coordinates.FaCentroids)';
FaCentroids = cell2mat(tmp);

%% main

% get components → take average per FA
[~, ~, listAveragesX] = intraFAanalysisExtract(maskCC, listIdx, targetX);
[~, ~, listAveragesY] = intraFAanalysisExtract(maskCC, listIdx, targetY);

orientationAvgTFM = [listAveragesX' listAveragesY'];

% calculate orientation wrt vector to center
% loop through FAs

tmp = cell2mat(vectorProx2Distal');

% ! note that both vectors life in different coordinate frames → flip! 
a = [];
a(:,1) = -tmp(:, 1); 
a(:,2) = -tmp(:, 2);
% a = flip(a, 2);
a(:,3) = 0;

b = [];
b(:,1) = orientationAvgTFM(:, 1);
b(:,2) = orientationAvgTFM(:, 2);
b = flip(b, 2);
b(:,3) = 0;

for i = 1:length(listIdx)
    % angleDeg(i) = atan2d(norm(cross(a(i,:),b(i,:))), dot(a(i,:),b(i,:)));
    % angleDeg(i) = atan2d( (  (a(i, 1) * b(i, 2)) - (a(i, 2) * b(i, 1))  ), (dot(a(i, 1:2), b(1,1:2))) );
    angleDeg(i) = acosd(dot(a(i,:) / norm(a(i,:)), b(i,:) / norm(b(i,:))));
    % angleDeg(i) = subspace(a(i,:), b(i,:));
end

% for i = 1:length(listIdx)
% %     angleDeg(i) = atan2d( (  (vectorProx2Distal{i}(2) * orientationAvgTFM(i, 1)) - (vectorProx2Distal{i}(1) * orientationAvgTFM(i, 2))  ), (dot(vectorProx2Distal{i}, orientationAvgTFM(1,1:2))) );
% angleDeg(i) = atan2d(norm(cross(a(i,:),b(i,:))), dot(a(i,:),b(i,:)));
% end

%% generate polarhistogram

nBins = 15; % 12: this generates 15 bins between +90° and -90° of 15° each

if showFigures
    figure
    h = polarhistogram(deg2rad(angleDeg), nBins);
else
end

%% validate using quiver plot

if showFigures

    % Reminder:
    % TFM array same axis as Matlab image
    % Matlab array →    Matlab image
    % x: vert down →    x: hor right
    % y: hor right →    y: vert down
    % **so simply flip x and y

    bw = ones(maskCC.ImageSize);
    bw = bw .* 220;
    bw = cast(bw, 'uint8');

    index = cell2mat(maskCC.PixelIdxList');
    bw(index) = 50;

    figure
    imshow(bw)
    hold on

    % validate TFM vectors
    X = FaCentroids(:, 1);
    Y = FaCentroids(:, 2);
    U = orientationAvgTFM(:, 1);
    V = orientationAvgTFM(:, 2);
    q = quiver(X,Y,U,V);
    q.Color = "b";
    q.LineWidth = 2;
    hold on

    % validate calcuated angle vectors
    sgn = -1;
    X = FaCentroids(:, 1);
    Y = FaCentroids(:, 2);
    U = sgn.*tmp(:, 2);
    V = sgn.*tmp(:, 1);
    q = quiver(X,Y,U,V);
    q.Color = "m";
    q.LineWidth = 2;
    hold on

    % angle validation
    angl = 45;
    ind = angleDeg > angl;
    p = plot(X(ind),Y(ind));
    p.LineStyle = "none";
    p.Marker = 'o';
    p.MarkerFaceColor = "k";
    p.MarkerEdgeColor = "none";
    p.MarkerSize = 7;
    hold on
    p = plot(X(ind),Y(ind));
    p.LineStyle = "none";
    p.Marker = 'o';
    p.MarkerFaceColor = "c";
    p.MarkerEdgeColor = "none";
    p.MarkerSize = 5;

    if sgn == -1
        orientation = 'Distal → Proximal vector';
    else
        orientation = 'Proximal → Distal vector';
    end

    title(['Blue: TFM tractions | Magenta: ' orientation ' | Cyan dots: |\theta| difference > ' num2str(angl) '°']);
else
end

%% log

% OneNote 2024-10-01 11:14 TFM vectors 

end


