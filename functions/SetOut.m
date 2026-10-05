function [outPath] = SetOut(outPath, inPath)
%SETOUT Summary of this function goes here
%   Detailed explanation goes here

outName = 'FALCON Output';

if (outPath == "")
    mkdir([inPath filesep outName])
    outPath = [inPath filesep outName];
else
end

