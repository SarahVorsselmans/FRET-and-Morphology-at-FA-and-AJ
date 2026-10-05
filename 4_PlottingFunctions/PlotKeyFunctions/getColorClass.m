function c = getColorClass(idx)
% RGB normalized to the range [0–1] (e.g. RGB of 71 becomes 71/255 = ~0.278)
    palette = [
        0.278, 0.376, 0.914;   % blue-purple-ish
        0.302, 0.686, 0.290;   % 2  green
        0.961, 0.796, 0.051;   % 7  yellow
        0.506, 0.278, 0.631;   % 4  purple
        0.855, 0.145, 0.114;   % 5  red
        0.102, 0.620, 0.667;   % 6  teal
        0.890, 0.467, 0.106;   % 3  orange
        0.647, 0.380, 0.204;   % 8  brown
        0.400, 0.400, 0.400;   % 9  grey
        0.851, 0.373, 0.608;   % 10 pink
        0.173, 0.447, 0.698;   % 1  blue
    ];
    c = palette(idx, :);
end

% function c = getColorClass(idx)
%     palette = [
%         0.173, 0.447, 0.698;   % 1  blue
%         0.302, 0.686, 0.290;   % 2  green
%         0.890, 0.467, 0.106;   % 3  orange
%         0.506, 0.278, 0.631;   % 4  purple
%         0.855, 0.145, 0.114;   % 5  red
%         0.102, 0.620, 0.667;   % 6  teal
%         0.961, 0.796, 0.051;   % 7  yellow
%         0.647, 0.380, 0.204;   % 8  brown
%         0.400, 0.400, 0.400;   % 9  grey
%         0.851, 0.373, 0.608;   % 10 pink
%     ];
%     c = palette(idx, :);
% end