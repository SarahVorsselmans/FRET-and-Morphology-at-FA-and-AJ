function saveVectorFigure(figHandle, filename, format, options)
% saveVectorFigure  Saves a figure as a vector file (PDF/SVG/EPS).
%
% Usage:
%   saveVectorFigure(gcf, 'output/myplot', 'pdf')
%   saveVectorFigure(gcf, 'plot.svg', 'svg', Width=12, Height=8)
%   saveVectorFigure(gcf, 'plot.pdf', 'pdf', Width=6.5, Height=4, Units='inches')
%
% Inputs:
%   figHandle    : Figure handle (e.g. gcf)
%   filename     : Output filename (with or without extension)
%   format       : 'pdf', 'svg', 'eps'
%
% Optional name-value inputs:
%   Width        : Figure width  (default: keep current)
%   Height       : Figure height (default: keep current)
%   Units        : 'centimeters' (default) | 'inches' | 'pixels' | 'points'
%
% Output:
%   Creates a vector graphic file on disk.

arguments
    figHandle           (1,1) matlab.ui.Figure
    filename            (1,:) char
    format              (1,:) char {mustBeMember(format, {'pdf','svg','eps'})}
    options.Width       (1,1) double = 0
    options.Height      (1,1) double = 0
    options.Units       (1,:) char  {mustBeMember(options.Units, ...
                            {'centimeters','inches','pixels','points'})} = 'centimeters'
end

% Strip extension if the user included one
[filepath, name, ~] = fileparts(filename);
outFile = fullfile(filepath, name + "." + format);

% Ensure parent folder exists
if ~isempty(filepath) && ~isfolder(filepath)
    mkdir(filepath);
end

% --- Resize figure if dimensions were provided ---
if options.Width > 0 || options.Height > 0
    prevUnits = figHandle.Units;
    figHandle.Units = options.Units;
    currentPos = figHandle.Position;

    newW = currentPos(3);
    newH = currentPos(4);
    if options.Width  > 0, newW = options.Width;  end
    if options.Height > 0, newH = options.Height; end

    figHandle.Position(3:4) = [newW, newH];
    figHandle.Units = prevUnits;
end

% Set figure renderer to painters for vector output
set(figHandle, 'Renderer', 'painters');

% Export using MATLAB's modern exportgraphics when possible
switch format
    case 'pdf'
        exportgraphics(figHandle, outFile, 'ContentType', 'vector');
    case 'svg'
        % SVG not supported by exportgraphics in some MATLAB versions
        print(figHandle, outFile, '-dsvg', '-vector');
    case 'eps'
        print(figHandle, outFile, '-depsc', '-vector');
end

fprintf('Saved vector figure: %s\n', outFile);
end

% function saveVectorFigure(figHandle, filename, format)
% % saveVectorFigure  Saves a figure as a vector file (PDF/SVG/EPS).
% %
% % Usage:
% %   saveVectorFigure(gcf, 'output/myplot', 'pdf')
% %   saveVectorFigure(gcf, 'plot.svg', 'svg')
% %
% % Inputs:
% %   figHandle : Figure handle (e.g. gcf)
% %   filename  : Output filename (with or without extension)
% %   format    : 'pdf', 'svg', 'eps'
% %
% % Output:
% %   Creates a vector graphic file on disk.
% 
%     arguments
%         figHandle  (1,1) matlab.ui.Figure
%         filename   (1,:) char
%         format     (1,:) char {mustBeMember(format, {'pdf','svg','eps'})}
%     end
% 
%     % Strip extension if the user included one
%     [filepath,name,~] = fileparts(filename);
%     outFile = fullfile(filepath, name + "." + format);
% 
%     % Ensure parent folder exists
%     if ~isempty(filepath) && ~isfolder(filepath)
%         mkdir(filepath);
%     end
% 
%     % Set figure renderer to painters for vector output
%     set(figHandle, 'Renderer', 'painters');
% 
%     % Export using MATLAB's modern exportgraphics when possible
%     switch format
%         case 'pdf'
%             exportgraphics(figHandle, outFile, 'ContentType','vector');
% 
%         case 'svg'
%             % Export using print (SVG not supported by exportgraphics in some MATLAB versions)
%             print(figHandle, outFile, '-dsvg', '-vector');
% 
%         case 'eps'
%             print(figHandle, outFile, '-depsc', '-vector');
%     end
% 
%     fprintf('Saved vector figure: %s\n', outFile);
% end