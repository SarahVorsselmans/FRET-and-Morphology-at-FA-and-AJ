function comparisonResults = Calc_compareConditionsToReference_ALL(CatTables, metrics, options)

arguments
    CatTables               struct
    metrics                 cell    = {}
    options.refDisease      string = 'CT'
    options.refExtra        string = 'rest'
end

    % Get all field names (e.g. 'TCatAll', 'TCat1', ...)
    tableNames = fieldnames(CatTables);

    % Initialize output struct
    comparisonResults = struct;

    % Loop over each table
    for i = 1:numel(tableNames)
        fname = tableNames{i};

        % Extract table
        inputTable = CatTables.(fname);

        % Run your existing function
        comparisonResults.(fname) = ...
            Calc_compareConditionsToReference(inputTable, metrics, ...
            refDisease = options.refDisease, ...
            refExtra   = options.refExtra);
    end
end