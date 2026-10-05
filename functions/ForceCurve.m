function [x_out, y_out, maxX, maxY] = ForceCurve(maxForce, sigma, maxFRET, deltaFRET, mu, showFlag, flagNorm)
%FORCECURVE Summary of this function goes here
%   Detailed explanation goes here

% Generate Function-value pairs
x = linspace(0,maxForce, 1000);
y = normpdf(x, mu, sigma); % generate gaussian distribution
y = y ./ max(y, [], 'all'); % normalize amplitude of distribution to [0 1]
y = y * (maxFRET - deltaFRET);
y = y + deltaFRET;

x_out = x;
y_out = y;

maxX = max(x, [], "all");
maxY = max(y, [], "all");

if flagNorm
    x = x / maxX;
    y = y / maxY;
else
end

if showFlag
    % Plot Gaussian
    figure("Name", 'Force Calibration curve')
    plot(x, y, "LineWidth", 3);
    hold on
    xlabel('Force (pN)', FontWeight='bold')
    ylabel('FRET Efficiency (%)', FontWeight='bold')
    title('Force Calibration VinTS7')
    hold on
else
end

end

