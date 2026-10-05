function [ExportT] = ExportTGraphPad(T, valueVar)
%EXPORTTGRAPHPAD Build ExportT from sorted T for GraphPad (minimal version)
%   ExportT = ExportTGraphPad(T, valueVar)
%   - T must contain: SampleFolder (group), EntryNumber (sort), and valueVar (e.g., 'FRETav')
%   - Groups follow current order in T ('stable'), within-group sort is ascending by EntryNumber.

    if nargin < 2 || isempty(valueVar)
        valueVar = 'FRETav'; % default to your current field
    end

    % ---- Basic validation (lightweight, but helpful) ----
    mustHave = {'SampleFolder', 'EntryNumber', valueVar};
    for i = 1:numel(mustHave)
        if ~ismember(mustHave{i}, T.Properties.VariableNames)
            error('ExportTGraphPad:MissingVar', 'Variable "%s" not found in T.', mustHave{i});
        end
    end
    if ~isnumeric(T.(valueVar))
        error('ExportTGraphPad:ValueNotNumeric', ...
              'The value variable "%s" must be numeric.', valueVar);
    end

    % ---- Get unique SampleFolders in current order ----
    sampleFolders = unique(string(T.SampleFolder), 'stable');
    nFolders = numel(sampleFolders);

    % ---- Collect valueVar per SampleFolder, sorted by EntryNumber ----
    valuesPerFolder = cell(nFolders,1);
    for k = 1:nFolders
        rows = string(T.SampleFolder) == sampleFolders(k);
        [~, idx] = sort(T.EntryNumber(rows), 'ascend');
        vals = T.(valueVar)(rows);
        valuesPerFolder{k} = vals(idx);
    end

    % ---- Determine max length and build wide matrix with NaN padding ----
    lenPerFolder = cellfun(@numel, valuesPerFolder);
    maxLen = max(lenPerFolder);
    W = nan(maxLen, nFolders);
    for j = 1:nFolders
        W(1:lenPerFolder(j), j) = valuesPerFolder{j};
    end

    % ---- Safe variable names for MATLAB ----
    safeNames = matlab.lang.makeValidName(cellstr(sampleFolders));

    % ---- Build ExportT table ----
    ExportT = array2table(W, 'VariableNames', safeNames);
    ExportT.EntryNumber = (1:maxLen).';
    ExportT = movevars(ExportT, 'EntryNumber', 'Before', 1);

    fprintf('ExportT created: %d rows x %d columns.\n', maxLen, nFolders);
end
