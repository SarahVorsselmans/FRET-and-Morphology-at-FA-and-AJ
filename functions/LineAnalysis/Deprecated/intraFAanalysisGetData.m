function [sourceFA, sourceFAmean] = intraFAanalysisExtract(lineFA, source)
%INTRAFAANALYSISEXTRACT Summary of this function goes here
%   Detailed explanation goes here

%     Copyright (C) 2024 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2024 Laurens Kimps
%

%% Get data from source 

% initialize variables
source = {};

% loop through FA segments
for i_FA = 1:length(lineFA.segmentPixels)
    for i_line = 1:length(lineFA(i_FA).segmentPixels)
        source(i_FA).values{i_line} = source(lineFA(i_FA).segmentPixels{i_line});
        source(i_FA).mean(i_line) = mean(source(i_FA).values{i_line});
    end
end

end

