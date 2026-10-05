function [outputList] = AveragePerRegion(intensity, connComp)
%AVERAGEPERREGION Summary of this function goes here

%     Copyright (C) 2024 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2024 Laurens Kimps
%     
%     <intensity>: The <connComp> mask will be applied to this image
%     <connComp>: a bwconncomp output 
%     <outputList>: the output vector/list containing the mean values per region

%% main
statsMeanIntensity = regionprops(connComp, intensity, 'MeanIntensity');
[averageIndexPerArea] = deal([statsMeanIntensity.MeanIntensity]);
outputList = averageIndexPerArea'; % rotate vector

end

