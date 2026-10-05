function [PaddedData] = LineAnalysisPadding(AveragedData, InterpolatedData, PaddingSize)
%LINEANALYSISPADDING Summary of this function goes here
%   Detailed explanation goes here


% Parameters
numFAs = numel(AveragedData);   % number of FAs
targetLength = PaddingSize;                       % desired number of columns

% Preallocate
PaddedData = zeros(numFAs, targetLength);

% Loop through FAs
for i = 1:numFAs
    lineData = AveragedData{i};      % 1 x N for this FA
    len = numel(lineData);

    if len <= targetLength
        % FA shorter than target length → pad with zeros
        PaddedData(i,1:len) = lineData;
        % remaining columns already 0 from preallocation
    else
        % FA longer than target length → use interpolated version
        % FA.intra.interp.FRET is a matrix (rows = FAs, columns = targetLength)
        PaddedData(i,:) = InterpolatedData(i,:);
    end
end
