function ShowFaIndex(faList)
%SHOWFAINDEX Summary of this function goes here
%   Detailed explanation goes here

%% parse data

mask = faList;

%% Main

imgLbl = zeros([faList.ImageSize]);
for i = 1:faList.NumObjects
    imgLbl(faList.PixelIdxList{i}) = i;
end

%% Draw figure

% make figure that tells you which FA is which with the data cursor
figure("Name", 'Index of FAs')
imshow(imgLbl)
title('If you hover over the FA you can see its index = FA label (use "data tips" tool)')

end

