function resultTable = Filter_IncludeExclude(allData, includePatterns, excludePatterns)
% Filter_ByPatterns
% -------------------------------------------------------------
% Filters rows of allData based on include and exclude patterns,
% applied to the SampleFolder column.
%
% INPUT:
%   allData          - table with column SampleFolder
%   includePatterns  - cell array of substrings to include
%   excludePatterns  - cell array of substrings to exclude
%
% OUTPUT:
%   resultTable      - filtered table
%
% Notes:
% - Matching is case‑sensitive (same as contains).
% - A row is included if:
%        (contains ANY include pattern) AND (contains NO exclude pattern)
% - If includePatterns is empty, ALL rows are included.
% - If excludePatterns is empty, no rows are excluded.
% -------------------------------------------------------------

    % Ensure SampleFolder is string array
    sampleCol = allData.SampleFolder;
    % if iscell(sampleCol)
    %     sampleCol = string(sampleCol);
    % elseif ischar(sampleCol)
    %     sampleCol = string(sampleCol);
    % end

    % ---- INCLUDE MASK ----
    if isempty(includePatterns)
        includeMask = true(height(allData), 1);   % include everything
    else
        includeMask = false(height(allData), 1);
        for i = 1:numel(includePatterns)
            includeMask = includeMask | contains(sampleCol, includePatterns{i});
        end
    end

    % ---- EXCLUDE MASK ----
    if isempty(excludePatterns)
        excludeMask = false(height(allData), 1);
    else
        excludeMask = false(height(allData), 1);
        for i = 1:numel(excludePatterns)
            excludeMask = excludeMask | contains(sampleCol, excludePatterns{i});
        end
    end

    % ---- Combine ----
    mask = includeMask & ~excludeMask;
    resultTable = allData(mask, :);
end

%% OLD CODE

% includePatterns = Info.includePatterns_Vin;
% excludePatterns = Info.excludePatterns_Vin;
% 
% includeMask = false(height(allData),1);
% for i = 1:numel(includePatterns)
%     includeMask = includeMask | contains(allData.SampleFolder, includePatterns{i});
% end
% 
% excludeMask = false(height(allData),1);
% for i = 1:numel(excludePatterns)
%     excludeMask = excludeMask | contains(allData.SampleFolder, excludePatterns{i});
% end
% 
% resultTable = allData(includeMask & ~excludeMask, :);
% T = resultTable;   % use filtered table (+ keep resultTable as backup)