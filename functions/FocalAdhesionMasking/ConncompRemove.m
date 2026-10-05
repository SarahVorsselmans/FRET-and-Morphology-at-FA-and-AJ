function [connComp] = ConncompRemove(connComp, delIdx)
%CONNCOMPREMOVE Summary of this function goes here

%     Copyright (C) 2024 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2024 Laurens Kimps
%     
%     <connComp>: bwconncomp output 

%% main

% remove PixelIdxList entries to be deleted
connComp.PixelIdxList(delIdx) = [];
% update the remaining number of objects
connComp.NumObjects = length(connComp.PixelIdxList);

end

