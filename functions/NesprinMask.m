function [smoothMask] = NesprinMask(img)
%NESPRINMASK Summary of this function goes here
%   Detailed explanation goes here
%% Parameters
out = 1;
in = 0;

strel_out = strel('disk',out);
strel_in = strel('disk',in);
%% LOG

%%% ----- Control

LOG.BW = edge(img,'log', 0, 6);
figure, imshow(LOG.BW)%, title('LOG - Binary Gradient Mask')


LOG.BWdfill = imfill(LOG.BW,'holes');
figure, imshow(LOG.BWdfill)%, title('LOG - Binary Image with Filled Holes')

[LOG.BWsingle] = SingleArea(LOG.BWdfill);
figure, imshow(LOG.BWsingle)%, title('LOG - Binary Image with Cleared')

se = strel('disk',1);
LOG.BWperimeter = LOG.BWsingle - imerode(LOG.BWsingle,se);
figure, imshow(LOG.BWperimeter)%, title('LOG - Control - PERIMETER 1')

if out == 0 & in == 0
    % DN
else
LOG.BWperimeter = imdilate(LOG.BWsingle,strel_out) - imerode(LOG.BWsingle,strel_in);
end

perim_size = in + out;
if perim_size < 2
    figure, imshow(LOG.BWperimeter)%, title('LOG - Control - PERIMETER')
    % Don't smooth
elseif perim_size <= 3
    [LOG.BWsmoothed] = smoothBinary(LOG.BWperimeter,3);
    figure, imshow(LOG.BWsmoothed)%, title('LOG - TS - Smoothed, kernel: 3')
elseif perim_size > 3
    [LOG.BWsmoothed] = smoothBinary(LOG.BWperimeter,5);
    figure, imshow(LOG.BWsmoothed)%, title('LOG - TS - Smoothed, kernel: 5 ')
end

% smoothMask = LOG.BWsmoothed;
smoothMask = LOG.BWperimeter;
%% FUNCTIONS

function [BW_single] = SingleArea(BW_many)
CC = bwconncomp(BW_many);
stats = regionprops(CC,'area');
[areas] = deal([stats.Area]);
[~,indx] = max(areas);
BW_single = false(size(BW_many));
BW_single(CC.PixelIdxList{indx}) = true;
end

function [BW_out] = smoothBinary(BW_in,wSize)
h = ones(wSize,wSize) ./ wSize^2;
result = imfilter(BW_in,h);
result(result > 0.5) = 1;
BW_out = false(size(BW_in));
BW_out(result == true) = true;
end

end


