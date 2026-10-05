function diffTbl = Calc_DifferenceTable(tableT, varNames, level1, level2, options)

arguments
    tableT table
    varNames string
    level1 string
    level2 string

    options.XCat string = "extrainfoCat"
    options.compareCat string = "diseaseCat"
    options.BioRepCol string = "expCat"

    options.includePatterns cell = {}
    options.excludePatterns cell = {}

    options.ExtraCols string = ["springCat" "proteinCat" "sampletypeCat"]
end

T = tableT;

% Filter
T = Filter_IncludeExclude(T, options.includePatterns, options.excludePatterns);

% Keep only requested levels
T = T(ismember(string(T.(options.XCat)), [level1 level2]), :);

% Ensure numeric variables
for i = 1:numel(varNames)
    T.(varNames(i)) = double(T.(varNames(i)));
end

%% Grouping variables
groupVars = [ ...
    options.compareCat ...
    options.BioRepCol ...
    options.XCat ...
    options.ExtraCols ];

groupVars = unique(groupVars,'stable');

%% Mean + SD of technical replicates
avgTbl = groupsummary( ...
    T,...
    cellstr(groupVars),...
    ["mean","std"],...
    cellstr(varNames));

grp    = avgTbl.(options.compareCat);
bioRep = avgTbl.(options.BioRepCol);
lvl    = string(avgTbl.(options.XCat));

groups = unique(grp,'stable');

%% -------------------------------------------------------------
% Collect output rows
%% -------------------------------------------------------------

rows = {};

for g = 1:numel(groups)
    idxG = grp == groups(g);
    bioRepsG = unique(bioRep(idxG));

    for b = 1:numel(bioRepsG)
        idx1 = idxG & bioRep == bioRepsG(b) & lvl == level1;
        idx2 = idxG & bioRep == bioRepsG(b) & lvl == level2;

        if ~(any(idx1) && any(idx2))
            continue
        end

        srcRow = find(idx1,1);
        
        row = avgTbl(srcRow,:); % metadata
        for v = 1:numel(varNames)
            row.(varNames(v) + "_SD") = NaN;
        end

        % loop
        for v = 1:numel(varNames)

            vn = varNames(v);
            meanCol = "mean_" + vn;
            stdCol  = "std_" + vn;

            mean1 = avgTbl{idx1, char(meanCol)};
            mean2 = avgTbl{idx2, char(meanCol)};

            sd1 = avgTbl{idx1, char(stdCol)};
            sd2 = avgTbl{idx2, char(stdCol)};

            sd1(isnan(sd1)) = 0;
            sd2(isnan(sd2)) = 0;

            % Difference
            row{1, char(meanCol)} = mean2 - mean1;

            % Propagated SD
            row.(vn + "_SD") = sqrt(sd1.^2 + sd2.^2);
        end
        rows{end+1,1} = row;
    end
end
if isempty(rows)
    diffTbl = avgTbl([],:);
    for v = 1:numel(varNames)
        diffTbl.(varNames(v) + "_SD") = [];
    end
    return
end
diffTbl = vertcat(rows{:});

%% -------------------------------------------------------------
% Rename mean_* columns back to original variable names
%% -------------------------------------------------------------

for v = 1:numel(varNames)
    vn = varNames(v);
    meanCol = "mean_" + vn;
    idx = strcmp(diffTbl.Properties.VariableNames, meanCol);
    if any(idx)
        diffTbl.Properties.VariableNames{idx} = char(vn);
    end
end

%% Remove std_* columns from groupsummary

varsToRemove = strings(0);

for v = 1:numel(varNames)
    varsToRemove(end+1) = "std_" + varNames(v);
end

varsToRemove = intersect( ...
    varsToRemove,...
    string(diffTbl.Properties.VariableNames));

diffTbl(:,cellstr(varsToRemove)) = [];

%% Replace category label

newLabel = categorical( ...
    "Diff_" + level2 + "v" + level1);

diffTbl.(options.XCat) = repmat( ...
    newLabel,...
    height(diffTbl),1);

end
