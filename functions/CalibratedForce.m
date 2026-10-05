function [forceCalibrated] = CalibratedForce(x, y, data, maxForce, sigma, maxFRET, deltaFRET, mu)
%CALIBRATEDFORCE SummargaussNorm of this function goes here
%   Detailed explanation goes here


for i = 1:length(data)
    [~,idx]=min(abs(y-data(i)));
    forceCalibrated(i) = x(idx);
end

