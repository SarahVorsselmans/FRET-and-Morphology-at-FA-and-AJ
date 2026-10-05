function tableT = Order_TableByCategories(tableT, orderingColumns, options)
% REORDERTABLEBYCATEGORIES Reorder table rows by unique combinations of category columns.
%
% INPUTS:
%   tableT           - Input table with category columns
%   orderingColumns  - Cell array of column names to combine (e.g. {'diseaseCat','extrainfoCat'})
%
% OPTIONAL INPUTS (name-value):
%   GroupOrder       - Cell array of group labels in desired order
%                      (e.g. {'CT_rest','CCM_rest','CCM_Piezo1',...})
%                      If omitted, groups are sorted alphabetically.
%
% OUTPUT:
%   tableT           - Table with rows reordered according to GroupOrder
%
% EXAMPLE:
%   T = reorderTableByCategories(T, {'diseaseCat','extrainfoCat'}, ...
%       'GroupOrder', {'CT_rest','CCM_rest','CCM_Piezo1'});

    arguments
        tableT          table
        orderingColumns (1,:) cell
        options.GroupOrder (1,:) cell = {}
    end

    % --- Build combined group label for each row ---
    % Concatenate the string values of each ordering column with '_'
    nRows   = height(tableT);
    nCols   = numel(orderingColumns);
    parts   = cell(nRows, nCols);

    for c = 1:nCols
        col = tableT.(orderingColumns{c});
        if iscategorical(col)
            parts(:, c) = cellstr(col);
        elseif iscell(col)
            parts(:, c) = col;
        elseif isstring(col)
            parts(:, c) = cellstr(col);
        else
            error('Column "%s" must be categorical, string, or cell array of chars.', ...
                  orderingColumns{c});
        end
    end

    % Join columns: 'CT' + 'rest'  ->  'CT_rest'
    rowLabels = join(string(parts), '_', 2);   % Nx1 string array
    rowLabels = cellstr(rowLabels);

    % --- Determine GroupOrder ---
    uniqueGroups = unique(rowLabels, 'stable');   % preserve first-occurrence order

    if isempty(options.GroupOrder)
        % Default: alphabetical
        groupOrder = sort(uniqueGroups);
    else
        groupOrder = options.GroupOrder;

        % Warn about groups in the table that are not listed in GroupOrder
        missing = setdiff(uniqueGroups, groupOrder);
        if ~isempty(missing)
            warning('reorderTableByCategories:missingGroups', ...
                'These groups exist in the table but are not in GroupOrder (appended at end):\n  %s', ...
                strjoin(missing, ', '));
            groupOrder = [groupOrder, missing(:)'];   % append unlisted groups at end
        end

        % Warn about groups in GroupOrder that don't exist in the table
        phantom = setdiff(groupOrder, uniqueGroups);
        if ~isempty(phantom)
            warning('reorderTableByCategories:phantomGroups', ...
                'These groups are in GroupOrder but not found in the table:\n  %s', ...
                strjoin(phantom, ', '));
        end
    end

    % --- Build the new row index ---
    newIdx = [];
    for g = 1:numel(groupOrder)
        mask   = strcmp(rowLabels, groupOrder{g});
        newIdx = [newIdx; find(mask)];  %#ok<AGROW>
    end

    % --- Apply reordering ---
    tableT = tableT(newIdx, :);
end