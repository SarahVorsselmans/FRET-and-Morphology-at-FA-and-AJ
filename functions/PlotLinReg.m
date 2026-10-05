function PlotLinReg(x, y, lineType, lineWidth, colorCode)
%PLOTLINREG Plot linear regression on current figure window
%   Detailed explanation goes here

    p = polyfit(x, y, 1);
    f = polyval(p ,x);
    [sortedX, srt] = sort(x);
    f = f(srt);
    plot(sortedX, f, lineType, "LineWidth", lineWidth, "Color", colorCode)
    hold on

end

