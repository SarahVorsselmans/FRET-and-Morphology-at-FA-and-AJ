%    This file is part of TFMcalc.
%     Copyright (C) 2013 Tecnun, School of Engineering, University of Navarra (Spain)
%     Copyright (C) 2016 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2012-2016 Alvaro Jorge-Penas
%
%     This library is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published
%     by the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This software is provided "as is",
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details
%     <http://www.gnu.org/licenses/>.


function [shiftBeads, corrFixedIm, corrMovingIm, corrCropThisImage] = StageShiftCorr_2(fixedIm, movingIm, shiftBeads, cropThisImage)
% This function corrects the misalignments between two bead images caused
% by innacurate stage shifts and other causes.
% The shift correction only affects the tranlational motion (i.e. no
% rotation, no higher order deformations)

fixedIm = dip_image(fixedIm);
movingIm = dip_image(movingIm);
cropThisImage(cropThisImage == 1) = 255;
cropThisImage = cast(cropThisImage, 'double');
cropThisImage = dip_image(cropThisImage);

sIm = size(movingIm);
beadDim = length(sIm);

shiftSmoothSigma = 3; % sigma for the gaussina smooting when computing the shift

% Change the boundary condtions to padd with zeros the locations where the
% image won't exist due to the shift correction
oldBoundaryOption = dip_getboundary(1);
dip_setboundary('add_zeros');

corrMovingIm = resample(movingIm,1,-shiftBeads,'bspline'); % the sign has to be inverted to be used with the "resample" function
% corrMovingIm = shift(movingIm,shiftBeads);

% Change the boundary conditions back to its previous value
dip_setboundary(oldBoundaryOption);

intShiftBeads = sign(shiftBeads).*ceil(abs(shiftBeads));

% Shift correction along the Z-axis
thOverlay = (1*thresholddip(fixedIm)) + (1*thresholddip(corrMovingIm));
if (beadDim==3)&&((100*(abs(intShiftBeads(3))/sIm(3)))>45)
    nonOverlayingCount =  zeros(1,sIm(3));
    for nn=0:sIm(3)-1
        slice = thOverlay(:,:,nn);
        nonOverlayingCount(nn+1) = sum(slice==1);
        clear slice
    end
    clear nn
    [~,cropIndx] = max(abs(diff(nonOverlayingCount)));
    nonOverlaying_ini = sum(thOverlay(:,:,0:cropIndx)==1);
    nonOverlaying_fini = sum(thOverlay(:,:,cropIndx+1:end)==1);
    if nonOverlaying_ini > nonOverlaying_fini % crop ini (top)
        if intShiftBeads(3)<0
            intShiftBeads(3) = sIm(3)+ intShiftBeads(3);
        end
    else % crop fini (bottom)
        if intShiftBeads(3)>0
            intShiftBeads(3) = -(sIm(3)-intShiftBeads(3));
        end
    end
    clear nonOverlayingCount cropIndx nonOverlaying_ini nonOverlaying_fini
end
clear thOverlay


% Crop the images accordingly
if (intShiftBeads(1)<=0)&&(intShiftBeads(2)<=0)
    if (beadDim==2)
        corrMovingIm = corrMovingIm(0:end+intShiftBeads(1),0:end+intShiftBeads(2));
        corrFixedIm = fixedIm(0:end+intShiftBeads(1),0:end+intShiftBeads(2));
        corrCropThisImage = cropThisImage(0:end+intShiftBeads(1),0:end+intShiftBeads(2));
        
    else
        if (intShiftBeads(3)<=0)
            corrMovingIm = corrMovingIm(0:end+intShiftBeads(1),0:end+intShiftBeads(2),0:end+intShiftBeads(3));
            corrFixedIm = fixedIm(0:end+intShiftBeads(1),0:end+intShiftBeads(2),0:end+intShiftBeads(3));
            corrCropThisImage = cropThisImage(0:end+intShiftBeads(1),0:end+intShiftBeads(2),0:end+intShiftBeads(3));
        else
            corrMovingIm = corrMovingIm(0:end+intShiftBeads(1),0:end+intShiftBeads(2),intShiftBeads(3):end);
            corrFixedIm = fixedIm(0:end+intShiftBeads(1),0:end+intShiftBeads(2),intShiftBeads(3):end);
            corrCropThisImage = cropThisImage(0:end+intShiftBeads(1),0:end+intShiftBeads(2),intShiftBeads(3):end);
        end
    end
    
    
elseif (intShiftBeads(1)<=0)&&(intShiftBeads(2)>0)
    if (beadDim==2)
        corrMovingIm = corrMovingIm(0:end+intShiftBeads(1),intShiftBeads(2):end);
        corrFixedIm = fixedIm(0:end+intShiftBeads(1),intShiftBeads(2):end);
        corrCropThisImage = cropThisImage(0:end+intShiftBeads(1),intShiftBeads(2):end);
    else
        if (intShiftBeads(3)<=0)
            corrMovingIm = corrMovingIm(0:end+intShiftBeads(1),intShiftBeads(2):end,0:end+intShiftBeads(3));
            corrFixedIm = fixedIm(0:end+intShiftBeads(1),intShiftBeads(2):end,0:end+intShiftBeads(3));
            corrCropThisImage = cropThisImage(0:end+intShiftBeads(1),intShiftBeads(2):end,0:end+intShiftBeads(3));
        else
            corrMovingIm = corrMovingIm(0:end+intShiftBeads(1),intShiftBeads(2):end,intShiftBeads(3):end);
            corrFixedIm = fixedIm(0:end+intShiftBeads(1),intShiftBeads(2):end,intShiftBeads(3):end);
            corrCropThisImage = cropThisImage(0:end+intShiftBeads(1),intShiftBeads(2):end,intShiftBeads(3):end);
        end
    end
    

    
elseif (intShiftBeads(1)>0)&&(intShiftBeads(2)<=0)
    if (beadDim==2)
        corrMovingIm = corrMovingIm(intShiftBeads(1):end,0:end+intShiftBeads(2));
        corrFixedIm = fixedIm(intShiftBeads(1):end,0:end+intShiftBeads(2));
        corrCropThisImage = cropThisImage(intShiftBeads(1):end,0:end+intShiftBeads(2));
    else
        if (intShiftBeads(3)<=0)
            corrMovingIm = corrMovingIm(intShiftBeads(1):end,0:end+intShiftBeads(2),0:end+intShiftBeads(3));
            corrFixedIm = fixedIm(intShiftBeads(1):end,0:end+intShiftBeads(2),0:end+intShiftBeads(3));
            corrCropThisImage = cropThisImage(intShiftBeads(1):end,0:end+intShiftBeads(2),0:end+intShiftBeads(3));
        else
            corrMovingIm = corrMovingIm(intShiftBeads(1):end,0:end+intShiftBeads(2),intShiftBeads(3):end);
            corrFixedIm = fixedIm(intShiftBeads(1):end,0:end+intShiftBeads(2),intShiftBeads(3):end);
            corrCropThisImage = cropThisImage(intShiftBeads(1):end,0:end+intShiftBeads(2),intShiftBeads(3):end);
        end
    end
    

    
elseif (intShiftBeads(1)>0)&&(intShiftBeads(2)>0)
    if (beadDim==2)
        corrMovingIm = corrMovingIm(intShiftBeads(1):end,intShiftBeads(2):end);
        corrFixedIm = fixedIm(intShiftBeads(1):end,intShiftBeads(2):end);
        corrCropThisImage = cropThisImage(intShiftBeads(1):end,intShiftBeads(2):end);
    else
        if (intShiftBeads(3)<=0)
            corrMovingIm = corrMovingIm(intShiftBeads(1):end,intShiftBeads(2):end,0:end+intShiftBeads(3));
            corrFixedIm = fixedIm(intShiftBeads(1):end,intShiftBeads(2):end,0:end+intShiftBeads(3));
            corrCropThisImage = cropThisImage(intShiftBeads(1):end,intShiftBeads(2):end,0:end+intShiftBeads(3));
        else
            corrMovingIm = corrMovingIm(intShiftBeads(1):end,intShiftBeads(2):end,intShiftBeads(3):end);
            corrFixedIm = fixedIm(intShiftBeads(1):end,intShiftBeads(2):end,intShiftBeads(3):end);
            corrCropThisImage = cropThisImage(intShiftBeads(1):end,intShiftBeads(2):end,intShiftBeads(3):end);  
        end
    end

    
end
        
% % Convert data to uint8 datatype
% corrFixedIm = uint8(round(stretch(corrFixedIm))); 
% corrMovingIm = uint8(round(stretch(corrMovingIm))); 

% Convert data to double datatype
corrFixedIm = double(corrFixedIm); 
corrMovingIm = double(corrMovingIm); 

corrCropThisImage = double(corrCropThisImage);


