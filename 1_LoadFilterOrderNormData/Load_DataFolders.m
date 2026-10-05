function allData = Load_DataFolders(Info, includePatterns)
% Load_DataFolders
% -------------------------------------------------------------
% Loops over multiple experiment folders, loads allOutputMetrics.mat
% inside each sample subfolder, extracts and cleans metrics,
% parses metadata from folder names, and returns a unified table.
%
% INPUTS:
%   mainDir          - Folder containing all experiment folders
%   dataFolders      - Cell array of experiment folder names
%   numFields        - Fields to convert to numeric
%   cellFields       - Fields to keep as cell arrays
%
% OUTPUT:
%   allData          - Combined table
%
% -------------------------------------------------------------

    mainDir     = Info.mainDir;
    dataFolders = Info.DataFolders;
    numFields   = Info.numFields; % for some reason saved as chr though
    cellFields  = Info.cellFields;
    
    % Initialize
    allTables = {};   % before loops
    idx = 1;

    for d = 1:numel(dataFolders)

        dataFolder = dataFolders{d};
        folderPath = fullfile(mainDir, dataFolder);

        % List sample folders
        sampleDirs = dir(folderPath);
        sampleDirs = sampleDirs([sampleDirs.isdir] & ~startsWith({sampleDirs.name}, '.'));

        % Filter sample folders (keep only those containing e.g. 'VinT')
        sampleDirs = sampleDirs(contains({sampleDirs.name}, includePatterns));

        for s = 1:numel(sampleDirs)

            sampleName = sampleDirs(s).name;
            matFile = fullfile(folderPath, sampleName, ...
                               'FALCON Output', 'allOutputMetrics.mat');

            if ~isfile(matFile)
                fprintf('Missing file: %s\n', matFile);
                continue;
            end

            fprintf('Loading %s...\n', matFile);
            S = load(matFile);

            if ~isfield(S, 'allOutputMetrics')
                warning('allOutputMetrics missing in %s', matFile);
                continue;
            end

            A = S.allOutputMetrics;
            nA = numel(A);
            if nA == 0
                continue;
            end

            %% ---------------------------------------------
            %   Extract metrics
            %% ---------------------------------------------
            metricData     = struct();
            metricDataAsIs = struct();
            validMask      = true(nA, 1);

            % ==== Numeric fields ==== make chr into number
            for m = 1:numel(numFields)
                f = numFields{m};

                vals = getFieldValues(A, f);
                if isempty(vals)
                    warning('%s not found in %s', f, matFile);
                    continue;
                end

                num = convertToNumeric(vals);
                validMask = validMask & ~isnan(num);
                metricData.(f) = num;
            end

            % ==== Cell (As‑is) fields ====
            for m = 1:numel(cellFields)
                f = cellFields{m};

                vals = getFieldValues(A, f);
                if isempty(vals)
                    warning('%s not found in %s', f, matFile);
                    continue;
                end
                metricDataAsIs.(f) = vals(:);
            end

            %% Apply mask
            numFieldNames = fieldnames(metricData);
            for f = numFieldNames'
                metricData.(f{1}) = metricData.(f{1})(validMask);
            end

            asIsFields = fieldnames(metricDataAsIs);
            for f = asIsFields'
                metricDataAsIs.(f{1}) = metricDataAsIs.(f{1})(validMask);
            end

            nEntries = sum(validMask);
            if nEntries == 0, continue; end

            %% ---------------------------------------------
            %   Build table for this sample
            %% ---------------------------------------------
            % Create base table first
            T = table(repmat({dataFolder}, nEntries, 1), ...
                      repmat({sampleName}, nEntries, 1), ...
                      (1:nEntries)', ...
                      'VariableNames', {'DataFolder','SampleFolder','EntryNumber'});

            % Add metadata columns by parsing
            T = addMetadataColumns(T, sampleName, nEntries);
            
            % Then append metrics
            T = [T, struct2table(metricData), struct2table(metricDataAsIs)];

            %% Store all rows
            allTables{idx} = T;
            idx = idx + 1;
            %allData = [allData; T];
        end
    end
    % Final concatenation
    if ~isempty(allTables)
        allData = vertcat(allTables{:});
    else
        allData = table();
    end
end


%% ============================================================
%    Helper: Extract a field from struct or cell array
%% ============================================================
function raw = getFieldValues(A, fieldName)

    raw = [];

    if isstruct(A) && isfield(A, fieldName)
        raw = {A.(fieldName)}';

    elseif iscell(A) && isstruct(A{1}) && isfield(A{1}, fieldName)
        raw = cellfun(@(x) x.(fieldName), A, 'UniformOutput', false);
        raw = raw(:);
    end
end

%% ============================================================
%    Helper: Convert values to numeric, or NaN if not possible
%% ============================================================
function num = convertToNumeric(raw)
    n = numel(raw);
    num = nan(n,1);

    for i = 1:n
        x = raw{i};

        if ischar(x) || (isstring(x) && isscalar(x))
            num(i) = str2double(x);
        elseif isnumeric(x) && isscalar(x)
            num(i) = x;
        else
            num(i) = NaN;
        end
    end
end

%% ============================================================
%    Helper: Parse metadata from sample folder name
%% ============================================================
function T = addMetadataColumns(T, sampleName, nEntries)

    tokens = regexp(sampleName, ...
        '(\d{8})_([^_]+)_([^_]+)_([^_]+)_?(.*)', ...
        'tokens', 'once');

    if isempty(tokens)
        T.date       = repmat({''}, nEntries, 1);
        T.sampletype = repmat({''}, nEntries, 1);
        T.spring     = repmat({''}, nEntries, 1);
        T.disease    = repmat({''}, nEntries, 1);
        T.extrainfo  = repmat({''}, nEntries, 1);
    else
        T.date       = repmat({tokens{1}}, nEntries, 1);
        T.sampletype = repmat({tokens{2}}, nEntries, 1);
        T.spring     = repmat({tokens{3}}, nEntries, 1);
        T.disease    = repmat({tokens{4}}, nEntries, 1);
        T.extrainfo  = repmat({tokens{5}}, nEntries, 1);
    end
end