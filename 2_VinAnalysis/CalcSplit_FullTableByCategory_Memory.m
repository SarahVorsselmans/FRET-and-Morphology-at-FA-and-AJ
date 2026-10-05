function TCat = CalcSplit_FullTableByCategory_Memory(T, Input, catIdx, options)
% CALCSPLIT_CATEGORYTABLE Processes a single category from T.
%
% USAGE:
%   TCat1 = CalcSplit_CategoryTable(T_Vin, Input, 1, FRET=true);
%   TCat2 = CalcSplit_CategoryTable(T_Vin, Input, 2, FRET=true);
%   TCat3 = CalcSplit_CategoryTable(T_Vin, Input, 3, FRET=true);
%   TCat4 = CalcSplit_CategoryTable(T_Vin, Input, 4, FRET=true);
%
% Then combine:
%   CatTables.TCat1 = TCat1;
%   CatTables.TCat2 = TCat2;
%   CatTables.TCat3 = TCat3;
%   CatTables.TCat4 = TCat4;

arguments
    T                   table
    Input               struct
    catIdx              (1,1) double
    options.FRET        logical = false
end

hasFRET = options.FRET;
nRows   = height(T);

Mask = cell(nRows, 1);
Int  = cell(nRows, 1);
if hasFRET
    FRETimageFA = cell(nRows, 1);
end

for rowIdx = 1:nRows
    mask         = double(T.MaskClass{rowIdx} == catIdx);
    Mask{rowIdx} = mask;

    tmp          = T.Int{rowIdx};
    tmp(mask==0) = NaN;
    Int{rowIdx}  = tmp;
    clear tmp

    if hasFRET
        tmp2                = T.FRETimageFA{rowIdx};
        tmp2(mask==0)       = NaN;
        FRETimageFA{rowIdx} = tmp2;
        clear tmp2
    end

    clear mask
end

% Assemble output table
TCat = table(T.DataFolder, T.SampleFolder, T.EntryNumber, ...
             T.diseaseCat, T.extrainfoCat, ...
             T.FA_summaryFrac(:, catIdx), ...
             Mask, Int, FRETimageFA, ...
             'VariableNames', ["DataFolder","SampleFolder","EntryNumber", ...
                               "diseaseCat","extrainfoCat","ClassFraction", ...
                               "Mask","Int","FRETimageFA"]);

% TCat = Calc_FA_AllMetrics(TCat, Input, FRET=hasFRET);
% TCat = Calc_FA_AveragingMetrics(TCat, Input, FRET=hasFRET);

end