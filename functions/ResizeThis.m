function [resized] = ResizeThis(input, sizeFinal)
%RESIZETHIS Summary of this function goes here
%   Detailed explanation goes here

logicalSize = size(input, 1);
sizeChange = round((sizeFinal - logicalSize) / 2);
resized = padarray(input,[sizeChange sizeChange], 0,'both');

end

