function ScatterColors(x,y)
%SCATTERCOLORS Summary of this function goes here
%   Detailed explanation goes here

colors = distinguishable_colors(length(x), [0 0 0] );

for i = 1:length(x)
scatter(x(i) , y(i) , 35,  'filled', "MarkerFaceColor", [colors(i,1), colors(i,2), colors(i,3) ]);
hold on



end

