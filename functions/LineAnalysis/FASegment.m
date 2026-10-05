function [bw_FA_mask] = FASegment(image, maskSingleCell, factor)  

imageSC = cast(image,'uint8');

imageSC = imageSC .* cast(maskSingleCell,'uint8');

% imageSC = wiener2(imageSC,[6 6]);
%   imageSC = adapthisteq(imageSC,'ClipLimit',0.01,'NumTiles',[15 15]);%lijkt niet voordelig?
% imageSC = imgaussfilt(imageSC,0.5);

threshold = graythresh(imageSC);
% threshold = threshold * 0.01;
imageBW = imbinarize(imageSC, threshold*factor);

% joinchannels('rgb', imageBW, imageSC)

% imageSC = imbinarize(imageSC,'adaptive','Sensitivity', 0.000001); % imageSC = imbinarize(imageSC,'adaptive','Sensitivity',0.0001);

CC = bwconncomp(imageBW);%making a list of all objects (FAs) in your FoV
L = labelmatrix(CC);%assign a different number to each object (FA) 
stats = regionprops(CC,'area');%calculating areas of each object (FA)
BW_ = ismember(L, find([stats.Area] > 2  )); %not taking the single pixels basically 
CC = bwconncomp(BW_);%new list of the objects 
%% next section removed the 'FAs' that are very low in intensity, probably they are not FA
statsMeanIntensity = regionprops(CC,image,'MeanIntensity');
[statsMeanIntensity] = deal([statsMeanIntensity.MeanIntensity]);
meanInt = mean([statsMeanIntensity],'all');
indx = find([statsMeanIntensity] > meanInt*0.5); %try putting 0 here 


numObj = length(indx);
% disp(numObj)
bw_FA_mask = false(size(imageBW));
for i = 1:numObj
    bw_FA_mask(CC.PixelIdxList{indx(i)}) = true;
end


end

