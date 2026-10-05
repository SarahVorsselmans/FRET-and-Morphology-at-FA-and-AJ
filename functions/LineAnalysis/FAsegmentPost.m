function [FASegmentoutput] = FAsegmentPost(Int, CellMask, minimumFA)
%FASEGMENTPOST Summary of this function goes here
%   Detailed explanation goes here


[bw_FA_mask1] = FASegment(Int, CellMask, 2.4);%takes the brighest
[bw_FA_mask2] = FASegment(Int, CellMask, 0.4);%takes the lowest
[bw_FA_mask3] = FASegment(Int, CellMask, 1.0);%takes in between

mask_rim = logical(bw_FA_mask2 - bw_FA_mask1);

skel = bwskel(mask_rim);

skel(bw_FA_mask1 == 1) = 0;

FAmask1 = bw_FA_mask3;
FAmask1(skel == 1) = 0;

CC = bwconncomp(FAmask1, 4);
L = labelmatrix(CC);
stats = regionprops(CC,'area');
BW = ismember(L, find([stats.Area] > 0  ));  %& [stats.Area] < medArea*60)
CC = bwconncomp(BW, 4);
L = labelmatrix(CC);
% dipshow(L)

% calculate outline/perimeter of FAs
outline = logical(imdilate(BW, strel('disk', 1)) - BW);
% dipshow(outline);

% filter based on size
CCoutline = bwconncomp(outline, 8);
Lout = labelmatrix(CCoutline);
stats = regionprops(CCoutline,'area');
% outline = ismember(Lout, find([stats.Area] > 4  )); % remove 'random' outline pixels
% dipshow(outline);
% include high threshold / exclude outline
BW2 = BW;
BW2(bw_FA_mask1 == 1) = 1;
BW2(outline == 1) = 0;

% close holes enclosed in FA area
se = strel('disk', 1); % Structural element for dilation/erosion
BW3 = imdilate(BW2, se); % Opening operation
BW3 = imerode(BW3, se); % Opening operation
BW3(outline == 1) = 0;
CC = bwconncomp(BW3, 4);
L3 = labelmatrix(CC);

% remove FAs that are too small
stats = regionprops(CC,'area');
BW3 = ismember(L3, find([stats.Area] > minimumFA  ));  %& [stats.Area] < medArea*60)

CC2 = bwconncomp(BW3, 4);
L2 = labelmatrix(CC2);
%   dipshow(L2); %use this if you want to adjust minimum FA threshold :)
finalValuestring = 'Adaptive';
% the idea: to remove skeleton when not needed
empty = zeros(CC2.ImageSize);
summedBW = empty;
for i = 1:CC2.NumObjects
    singleBW = empty;
    singleBW(cell2mat(CC2.PixelIdxList(i))) = 1;
    se = strel('disk', 1); % Structural element for dilation/erosion
    singleBW = imdilate(singleBW, se); % Opening operation
    summedBW = summedBW + singleBW; % overlapping areas will have a value greater than 1
end
skel2 = empty;
skel2(summedBW > 1) = 1;
newBW = bw_FA_mask3 - skel2;
newBW(outline == 1) = 0;
% filter based on size
CCnew = bwconncomp(newBW, 4);
Lnew = labelmatrix(CCnew);
stats = regionprops(CCnew,'area');
newBW = ismember(Lnew, find([stats.Area] > minimumFA )); % remove small FAs
% smallBW = ismember(Lnew, find([stats.Area] == 1 )); % remove small FAs
% newBW(smallBW == 1) = 0;
CCnew2 = bwconncomp(newBW, 4);
Lnew2 = labelmatrix(CCnew2);
dipshow(Lnew2)
% dipshow(skel2);
FASegmentoutput = newBW; % output this to the rest of the code


end

