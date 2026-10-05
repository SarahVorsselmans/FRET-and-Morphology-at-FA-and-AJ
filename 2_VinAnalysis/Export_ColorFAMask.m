function Export_ColorFAMask(T, options)
% Export_ColorFAMask
% --------------------------------------------------------
% Converts T.MaskClass to a colored RGB image and saves it
% as a .tif file for every row in T.
%
% INPUTS:
%   T       = main table with columns:
%               - MaskClass   : cell array of 2D double (values 0–N)
%               - Path        : char array with output folder per row
%               - EntryNumber : numeric ID per row
%   options = (reserved for future use)
%
% OUTPUT filename format: Cell<EntryNumber>_Class.tif
% --------------------------------------------------------

arguments
    T          table
    options.includePatterns cell      = {}
    options.excludePatterns cell      = {}
end

nCats = size(T.FA_summaryFrac, 2);  % calculates number of columns (= categories)
% check to see if nCats does not exceed number of colors
assert(nCats <= 10, ...
    'Export_ColorFAMask: nCats=%d exceeds the 10-color palette in getColorClass.', nCats);

    %% Filter rows via include/exclude on SampleFolder
    T = Filter_IncludeExclude(T, options.includePatterns, options.excludePatterns);

    %% --- Build colormap (5 x 3, for values 0–4) ---
    cmap = [0, 0, 0];  % background
    for c = 1:nCats
        cmap = [cmap; getColorClass(c)];
    end

    nRows = height(T);

    for i = 1:nRows

        %% --- Check path exists ---
        savePath = T.Path{i};
        if ~isfolder(savePath)
            warning('Export_ColorFAMask: folder does not exist for row %d, skipping.\n  Path: %s', ...
                i, savePath);
            continue;
        end

        %% --- Get classification mask ---
        ClassMask = T.MaskClass{i};
        if isempty(ClassMask)
            warning('Export_ColorFAMask: T.MaskClass{%d} is empty, skipping.', i);
            continue;
        end

        %% --- Convert indexed mask to RGB image ---
        % ind2rgb expects 1-based indices, so shift values from (0–4) to (1–5)
        ClassMaskIndexed = ClassMask + 1;
        RGBimage = ind2rgb(ClassMaskIndexed, cmap);  % returns H x W x 3 double [0,1]

        %% --- Build filename and full save path ---
        fileName  = sprintf('Cell%d_Class.tif', T.EntryNumber(i));
        fullPath  = fullfile(savePath, fileName);

        %% --- Save as TIFF ---
        imwrite(RGBimage, fullPath, 'tif');

    end
    fprintf('Done.\n');
end