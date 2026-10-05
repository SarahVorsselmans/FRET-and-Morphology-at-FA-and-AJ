function [cellArray, arrayAvgNoNAN, listAverages] = intraFAanalysisExtract(maskCC, listIdx, target)
%INTRAFAANALYSISEXTRACT Extracts per-line and per-FA averages
%   Leaves cellArray unchanged.
%   Removes resulting NaNs only from per-line averages for interpolation.

numFAs = length(listIdx);

cellArray     = cell(1, numFAs);  % store per-line raw pixel arrays
arrayAvgNoNAN = cell(1, numFAs);  % store per-FA averages with NaNs removed
listAverages  = NaN(1, numFAs);   % FA-level mean

for i = 1:numFAs
    allPixels = [];    % accumulate valid pixels for FA-level mean
    lineAverages = zeros(1, length(listIdx{i})); % preallocate

    % loop through each perpendicular line of this FA
    for k = 1:length(listIdx{i})
        % extract pixel values from target
        pixelVals = target(listIdx{i}{k});
        cellArray{i}{k} = pixelVals; % raw data stays untouched

        % calculate average of this line, ignoring NaNs
        lineAverages(k) = mean(pixelVals, "omitnan");

        % accumulate all valid pixels
        allPixels = [allPixels, pixelVals(:)']; %#ok<AGROW>
    end

    % remove NaNs from the resulting line averages
    arrayAvgNoNAN{i} = lineAverages(~isnan(lineAverages));

    % FA-level mean using all valid pixels
    validPixels = allPixels(~isnan(allPixels));
    if ~isempty(validPixels)
        listAverages(i) = mean(validPixels, "omitnan");
    else
        listAverages(i) = NaN; % FA had no usable pixels
    end
end
