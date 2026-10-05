function [s, sMean, syavg, syavgMean] = interpolateLineAnalysis(titleName, lineData, lineLength, user)
%INTERPOLATELINEANALYSIS Summary of this function goes here
%   Detailed explanation goes here

%% parse

outPath = user.outputPath;
saveName = user.saveName; 
%% init
test = 1;
addpath('functions\community\inpaint_nans\');

%% main

% get xq = x-axis values starting from 0 w/ length = <lineLength>
xq = linspace(0, lineLength-1, lineLength);
% normalize [0 1]
xq = xq ./ max(xq);
% we need xq later on for <s>

disp(newline)

% loop through FAs
for i = 1:length(lineData)

    % get x axis values starting from 0 w/ length = length(lineData{i})
    lengthLine = length(lineData{i});
    x = linspace(0, lengthLine-1, lengthLine);
    x = x ./ max(x);

    % get y axis data from line analysis
    y = lineData{i};
    % Normalize y
    yNorm = y ./ max(y);

    % get s, the interpolated version of y
    if nnz(isnan(y)) < 1 & length(y) >= 10
        y = inpaint_nans(y, 1); % see website for description
        disp(['     interpolating line data: #' num2str(i)])
        % get interpolation of <y>
        s(i, :) = interp1(x,y,xq, "linear");
        % get interpolation of averaged <y>
        syavg(i, :) = interp1(x, yNorm, xq, "linear");
    else
        s(i, :) = repmat(NaN, 1, lineLength);
        syavg(i, :) = repmat(NaN, 1, lineLength);
    end
end

sMean = mean(s, 1, "omitnan");
syavgMean = mean(syavg, 1,"omitnan");
%% for testing only

if 0
    s = s(~isnan(s));
end

%% test interp1

color1 = [0.4940 0.1840 0.5560];

if 0
    figure
    scatter(x, y)
    hold on
    scatter(xq, s(i))
else
end

if test
    figure
    title([titleName ' • Line analysis' ' (#FAs = ' num2str(length(lineData)) ')'])
    xlabel('← Proximal | FA | Distal →')
    ylabel('NOT normalized data')
    hold on
    for k = 1:length(lineData)
        p = plot(xq, s(k, :));
        hold on
        p.LineStyle = "-";
        p.Color = color1.* 1.6;
        p.LineWidth = 1;
    end
    hold on
    p2 = plot(xq, sMean);
    p2.LineWidth = 5;
    p2.Color = color1;
else
end

saveas(gcf,[outPath filesep saveName 'Line_' titleName '.png'])


if test
    figure
    title([titleName ' • Line analysis' ' (#FAs = ' num2str(length(lineData)) ')'])
    xlabel('← Proximal | FA | Distal →')
    ylabel('Normalized data')
    hold on
    for k = 1:length(lineData)
        p3 = plot(xq, syavg(k, :));
        p3.LineStyle = "-";
        p3.Color =  color1.* 1.6;
        hold on
        p3.LineWidth = 1;
    end
    hold on
    p4 = plot(xq, syavgMean);
    p4.LineWidth = 3;
    p4.Color = color1;
else
end

saveas(gcf,[outPath filesep saveName 'Line_' titleName '_Norm.png'])


end

