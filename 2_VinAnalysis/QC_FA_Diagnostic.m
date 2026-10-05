function QC_FA_Diagnostic(FRETmask, EccentricityValues, varargin)
% Plot_FA_Diagnostic
% -----------------------------------------------------------
% Creates a diagnostic overlay plot showing:
%   - original FRETindex mask (grayscale)
%   - detected FA boundaries (red)
%   - FA labels in green:
%          "<ID> (Ecc=<value>)"
%
% INPUTS:
%   FRETmask           : numeric 2D mask (background=100, FA ≠100)
%   EccentricityValues : Nx1 vector of eccentricity values
%                        (corresponds to stats.Eccentricity)
%
% OPTIONAL NAME-VALUE PAIRS:
%   'ShowLabels'       : true/false (default = true)
%
% -----------------------------------------------------------

p = inputParser;
p.addParameter('ShowLabels', true, @(x)islogical(x));
p.parse(varargin{:});
opt = p.Results;

% Defensive copy
M = FRETmask;

if isempty(M) || ~isnumeric(M)
    error('FRETmask must be a numeric 2D matrix.');
end

% Replace invalid values with background=100
M(~isfinite(M)) = 100;

% Binary FA mask
BW = (M ~= 100);

% Connected components
CC = bwconncomp(BW, 8);
stats = regionprops(CC, 'Centroid');

% Extract centroids just for positioning labels
centroids = vertcat(stats.Centroid);

% Prepare outline
outline = bwperim(BW);

figure; hold on;

% --- Show original mask ---
imagesc(M);
colormap(gray);
axis image off;
title('FA Segmentation Diagnostic');

% --- Overlay outlines ---
hOutline = imshow(outline);
set(hOutline, 'AlphaData', outline * 0.8);
set(hOutline, 'CData', cat(3, ones(size(outline)), zeros(size(outline)), zeros(size(outline))));  % red

% --- FA labels in green WITH eccentricity ---
if opt.ShowLabels && ~isempty(centroids)
    
    % EccentricityValues should match #FAs
    nFA = size(centroids,1);
    if nargin < 2 || isempty(EccentricityValues)
        EccentricityValues = nan(nFA,1);  % placeholder if missing
    end
    
    for k = 1:nFA
        eccStr = sprintf('%.2f', EccentricityValues(k));  % rounded
        label  = sprintf('%d (%s)', k, eccStr);        % label format
        
        text(centroids(k,1)+3, centroids(k,2), label, ...
            'Color','g','FontSize',8,'FontWeight','bold');
    end
end

hold off;

end