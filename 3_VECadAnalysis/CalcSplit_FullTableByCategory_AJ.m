function CatTables = CalcSplit_FullTableByCategory_AJ(T, Input, options)
% SPLITTABLEBYCATEGORY Splits table T into per-category tables based on T.MaskClass,
% then applies Calc_FA_AllMetrics and Calc_FA_AveragingMetrics to each category table.
%
% For each category index N (derived from Input.CatCategories):
%   - CatTables.TCatN.Mask        : binary mask (1 where MaskClass==N, 0 elsewhere)
%   - CatTables.TCatN.Int         : T.Int with non-object pixels set to NaN
%   - CatTables.TCatN.FRETimageFA : T.FRETimageFA with non-object pixels set to NaN
%                                   (only if T.FRETimageFA exists)
%
% INPUTS:
%   T       - Table with columns: MaskClass, Int, and optionally FRETimageFA.
%             Each cell contains a 512x512 double image.
%   Input   - Structure with field CatCategories (cell array, one entry per category).
%
% OUTPUT:
%   CatTables - Struct of tables, one per category: CatTables.TCat1, CatTables.TCat2, ...
%               Each table has metrics added by Calc_FA_AllMetrics and
%               Calc_FA_AveragingMetrics.

arguments
    T                       table
    Input                   struct
    options.FRET            logical = false
end

hasFRET      = options.FRET;
nCats        = numel(Input.AJCatCategories);
nRows        = height(T);

% To save memory
allmaskClass     = T.MaskClass;
allInt           = T.Int;
if hasFRET
    allFRET      = T.FRETimageFA;
end

    % -------------------------------------------------------------------------
    % Step 1: Build per-category tables
    % -------------------------------------------------------------------------
    for catIdx = 1:nCats

        % Pre-allocate cell columns for this category
        Mask = cell(nRows, 1);
        Int  = cell(nRows, 1);
        if hasFRET
            FRETimageFA = cell(nRows, 1);
        end
        % --- Transfer columns ---
        DataFolder     = T.DataFolder;
        SampleFolder   = T.SampleFolder;
        EntryNumber    = T.EntryNumber;
        %Path           = T.Path;
        %DateCat        = T.DateCat;
        diseaseCat     = T.diseaseCat;
        extrainfoCat   = T.extrainfoCat;
        ClassFraction  = T.FA_summaryFrac(:, catIdx);

        for rowIdx = 1:nRows
            % --- Build binary mask -------------------------------------------
            maskClass    = allmaskClass{rowIdx};          % 512x512 double
            binaryMask   = double(maskClass == catIdx);  % 1 where catIdx, 0 elsewhere
            Mask{rowIdx} = binaryMask;

            % --- Apply mask to Int -------------------------------------------
            nanMask               = binaryMask;
            nanMask(nanMask == 0) = NaN;                 % 0->NaN, 1 stays 1
            % nanMask(nanMask == 0) = 0;                 % Alt: keep as 0 instead of NaN

            Int{rowIdx} = allInt{rowIdx} .* nanMask;

            % --- Apply mask to FRETimageFA (if present) ----------------------
            if hasFRET
                FRETimageFA{rowIdx} = allFRET{rowIdx} .* nanMask;
            end

        end % rows

        % --- Assemble per-category table -------------------------------------
        if hasFRET
            CatTable = table(Mask, Int, FRETimageFA, DataFolder, SampleFolder, ...
                EntryNumber, diseaseCat, extrainfoCat, ClassFraction); % Path, DateCat
        else
            CatTable = table(Mask, Int, DataFolder, SampleFolder, ...
                EntryNumber, diseaseCat, extrainfoCat, ClassFraction); % Path, DateCat
        end

        fieldName             = sprintf('TCat%d', catIdx);
        CatTables.(fieldName) = CatTable;

        % Empty columns so it doesn't crash with CatAll
        CatTables.(fieldName).CellAvUmSize       = nan(nRows, 1);
        CatTables.(fieldName).CellAv_MajAx_um    = nan(nRows, 1);
        CatTables.(fieldName).CellAv_MinMajRatio = nan(nRows, 1);
        CatTables.(fieldName).CellAlignIndex     = nan(nRows, 1);

    end % categories

    % -------------------------------------------------------------------------
    % Step 2: Apply metric functions to each per-category table
    % -------------------------------------------------------------------------
    for catIdx = 1:nCats

        fieldName = sprintf('TCat%d', catIdx);

        CatTables.(fieldName) = Calc_Thickness(CatTables.(fieldName), Input);
        %CatTables.(fieldName) = Calc_FA_AveragingMetrics(CatTables.(fieldName), Input);

    end % categories

    %% Add original table to structure as CatAll
    CatTables.CatAll = T;
    CatTables.CatAll.ClassFraction = nan(nRows, 1); % to make sure next functions don't crash on this column

end
