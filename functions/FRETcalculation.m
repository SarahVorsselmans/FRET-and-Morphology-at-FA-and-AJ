function [FRET_image, FRETav,  FRET_list, Title, Plow, Phigh] = FRETcalculation(FRET_image, mask, user)
%FRETCALCULATION Summary of this function goes here

FRET_Threshold = user.FRET_Threshold;

FRET_image(mask == 0 | FRET_image <= 0 | FRET_image >= FRET_Threshold) = NaN;

FRET_list = FRET_image(~isnan(FRET_image) & FRET_image ~= 0);
Plow = round(prctile(FRET_list,1)); %prctile function creates dynamic numbers based on the intensity distribution
Phigh = round(prctile(FRET_list,99));

FRETav = num2str(mean(FRET_list));
Title = ['Average FRET efficiency = ' FRETav '%'];

end

