function longTbl = Help_expandMatrixColumn(tbl, matColName, levelNames)
% Converts an N×4 matrix column into long format, adding a 'fractionCat' column.
%
%   tbl        : original table
%   matColName : name of the N×4 double column, e.g. 'FA_summaryFrac'
%   levelNames : 1×4 cell of labels, e.g. {'bound','free','partial','other'}

    mat  = tbl.(matColName);           % N×4 double
    nCol = size(mat, 2);
    assert(nCol == numel(levelNames), 'levelNames length must match number of sub-columns');

    parts = cell(nCol, 1);
    for c = 1:nCol
        t = tbl;
        t.(matColName)  = mat(:, c);           % replace matrix with scalar column
        t.fractionCat   = repmat(string(levelNames{c}), height(t), 1);
        parts{c} = t;
    end
    longTbl = vertcat(parts{:});
end