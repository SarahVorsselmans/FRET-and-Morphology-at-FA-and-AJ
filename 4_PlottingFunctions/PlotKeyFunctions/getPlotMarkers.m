function markers = getPlotMarkers(n)
%GETPLOTMARKERS  Return a 1-by-n cell array of marker style strings.
%
%   markers = getPlotMarkers(n)
%
% Input
%   n       : number of marker styles needed (cycled if n > 12)
%
% Output
%   markers : 1-by-n cell array of marker character strings

if nargin < 1 || isempty(n), n = 1; end

markerSet = {'o','s','^','d','v','>','<','p','h','x','+','*'};
nM        = numel(markerSet);
markers   = markerSet(mod(0:n-1, nM) + 1);
end
