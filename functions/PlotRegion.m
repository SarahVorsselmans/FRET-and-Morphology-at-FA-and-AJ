function PlotRegion(cc, numObj, idxList, imgSize, figTitle, flagTesting, cnt, imgcnt)
%PLOTFA Summary of this function goes here
%   Detailed explanation goes here

% Generate color list for displaying
colors = distinguishable_colors(numObj, [0 0 0] );

% generate zero images
img = zeros([imgSize 3]);
img1 = zeros([imgSize]);
img2 = zeros([imgSize]);
img3 = zeros([imgSize]);
imgFull = zeros([imgSize]);
imgLbl = zeros([imgSize]);

% generate images (color and binary and label images)
for i = 1:numObj
    img1(idxList{i}) = colors(i,1); % r
    img2(idxList{i}) = colors(i,2); % g
    img3(idxList{i}) = colors(i,3); % b
    labelList{i} = num2str(i); % for the legend
    imgFull(idxList{i}) = true; % for the testing option
    imgLbl(idxList{i}) = i; % for the label picker figure
end
img(:,:,1:3) = cat(3, img1, img2, img3); % concatenate to create rgb image
subplot(1,1,imgcnt)
imshow(img)
text = ([num2str(size(img, 1)) ' x ' num2str(size(img, 2))]);
listLength = num2str(numObj);
title([figTitle listLength ', Resolution: ' text ')'])
% Add legend that tells you which FA has which unique number/label
hold on
for i = 1:numObj
    scatter([],[],10,[colors(i,1) colors(i,2) colors(i,3)],'filled');
end
legend(labelList)

% Check that the legend is correct
if flagTesting == 1 & cnt == 0

    i = ceil(numObj/2); % which FA do you want to check?
    imgTest = zeros([imgSize 3]);
    imgTest(idxList{i}) = true;
    imFused = imfuse(imgFull, imgTest);
    figure
    imshow(imFused)
    title(['This should be FA number ' num2str(i) ' (in the updated list)'])
    cnt = 1;

end

end

