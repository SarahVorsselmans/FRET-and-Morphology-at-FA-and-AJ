function [movingOut] = ShiftCorrChannel(fixed, moving)
%SHIFTCHANNEL Summary of this function goes here

%     Copyright (C) 2024 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2024 Laurens Kimps
%
%     This function makes use of stageShiftCorr 
%     Copyright (C) 2013 Tecnun, School of Engineering, University of Navarra (Spain)
%     Copyright (C) 2016 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2012-2016 Alvaro Jorge-Penas
%     
%     <fixed>: image to which <moving> is moved (square image)
%     <moving>: image which is moved and cropped and zero padded to match size of <fixed> again 
%
%     → if no shift: moving is directly put as unaltered output variable

%% main
% check spatial alignation of data
[shift, corrFixedIm, corrMovingIm] = stageShiftCorr(fixed, moving);

%     % Display shift in overlayed images (for testing/diagnosis)
    joinchannels('rgb', corrFixedIm, corrMovingIm) % corrected images !! with intensity rescaling !!
    joinchannels('rgb', fixed, moving) % uncorrected images

% if shift != [0;0] → add zero padding to correct for cropping
if any(shift)

    % move parameters
    y = abs(shift(1)); % flip x y (image vs array convention)
    x = abs(shift(2));
    sz = size(fixed, 1); % assume square images

    disp(['• Acceptor analysis: shift correcting acceptor channel (shift[x y] = ' num2str(shift') ')'])

    % create target array
    emptytarget = false(size(fixed));
    emptytarget = cast(emptytarget, class(fixed));

    % shift arrays
    if shift(1) && shift(2) >= 0
        emptytarget(1+x:sz, 1+y:sz) = moving(1:sz-x, 1:sz-y);
    elseif shift(1) && shift(2) <= 0
        emptytarget(1:sz-x, 1:sz-y) = moving(1+x:sz, 1+y:sz);
    elseif shift(1) >= 0 && shift(2) <= 0
        emptytarget(1+x:sz, 1:sz-y) = moving(1:sz-x, 1+y:sz);
    elseif shift(1) <= 0 && shift(2) >= 0
        emptytarget(1:sz-x, 1+y:sz) = moving(1+x:sz, 1:sz-y);
    else
    end

    %     % Testing/diagnosis of shift correction
    %     joinchannels('rgb', fixed, emptytarget) % corrected images !! no intensity rescaling !!

    % parse output
    movingOut = emptytarget;

else
    disp(['• Acceptor analysis: No shift correction needed (shift[x y] = ' num2str(shift') ')'])
    movingOut = moving;
end



end

