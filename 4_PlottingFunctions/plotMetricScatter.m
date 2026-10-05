function [fig, pVals, overviewTbl, anovaTbl] = plotMetricScatter(tableT, varName, options)
%PLOTMETRICSCATTER  Scatter plot of a numeric metric per x-category, split by
%                   one or more comparison groups. Points are coloured by
%                   DataFolder; group means are drawn as short horizontal lines;
%                   a gray vertical band highlights the last compareGroup side of
%                   each category tick; violin plots show the distribution behind
%                   the scatter.
%
%                   Statistical testing is performed on the unfiltered dataset
%                   so that p-values remain stable across figures.
%
% USAGE
%   [fig, pVals, overviewTbl, anovaTbl] = plotMetricScatter( ...
%       tableT, varName, Name, Value, ...)
%
% REQUIRED INPUTS
%   tableT    : MATLAB table containing all data columns.
%   varName   : char | string | numeric vector.
%               Name of the numeric metric column, or a pre-extracted numeric
%               vector (requires 'VarColName' + 'VarColIndex' as well).
%
% OPTIONAL NAME-VALUE INPUTS
%  --- Data selection ---
%   'levels'           : Cell array. Ordered x-axis levels from XCat.
%                        Auto-inferred from XCat column when empty.
%   'XCat'             : Column name for x-axis categories. Default: 'extrainfoCat'.
%   'compareCat'       : Column name that splits each x-level into groups.
%                        Default: 'diseaseCat'.
%   'compareGroups'    : Cell array. Subset/order of compareCat values to plot.
%                        Default: all unique non-empty values.
%   'ConditionOrder'   : Cell array. Explicit left-to-right tick ordering.
%                        Entries not present in data are silently dropped.
%   'includePatterns'  : Cell array. Keep rows whose SampleFolder matches any pattern.
%   'excludePatterns'  : Cell array. Drop rows whose SampleFolder matches any pattern.
%   'SplitCat'         : Column name for an additional pre-filter.
%   'SplitCatLevel'    : Value to retain in SplitCat.
%   'CutoffBottom'     : Numeric scalar. Values below are set to NaN.
%   'CutoffTop'        : Numeric scalar. Values above are set to NaN.
%
%  --- Sample-type expansion ---
%   'IncludeTL'        : Logical. When true, the 'sampletypeCat' column is
%                        incorporated so that each condition tick is labelled
%                        compareGroup + "_" + sampleType + "_" + level
%                        (e.g. CT_TS_rest, CT_TL_rest, CCM_TS_rest).
%                        TL entries are only present for levels where TL data
%                        actually exists; absent combos are dropped automatically.
%                        Default: true.
%   'SampleTypeCol'    : Name of the sample-type column. Default: 'sampletypeCat'.
%   'SampleTypeRef'    : The reference / always-present sample type (placed
%                        first within each level). Default: 'TS'.
%
%  --- Statistics ---
%   'statTest2'        : Primary statistical test.
%                        'lmm' (default) | 'nested' | 'onewayTech' |
%                        'onewayBio' | 'rm' | 'ttest'.
%   'Contrasts'        : Contrasts passed to the stat function.
%                        'all' (default) or explicit cell array of pairs.
%   'ContrastRemoveG1' : Cell array. Remove these labels from Group-1 side
%                        when building auto contrasts.
%   'ContrastRemoveG2' : Cell array. Remove these labels from Group-2 side.
%   'CollapseRep'      : Logical. Average technical replicates per bio-rep
%                        before testing. Default: true.
%   'TTestType'        : 'paired' (default) | 'unpaired'. Used when statTest2='ttest'.
%
%  --- Appearance ---
%   'grayband'         : Logical. Draw gray band behind last compareGroup. Default: true.
%   'Palette'          : Color palette name passed to getPlotColors. Default: 'distinct'.
%   'MultiplePalette'  : Logical. Use a different palette per condition combo. Default: true.
%   'condColors'
%   'MarkerSize'       : Scalar marker area. Default: 36.
%   'Alpha'            : Marker transparency [0,1]. Default: 0.75.
%   'Jitter'           : Horizontal jitter amplitude. Default: 0.08.
%   'ViolinWidth'      : Max half-width of each violin. Default: 0.13.
%   'ViolinAlpha'      : Violin transparency [0,1]. Default: 0.25.
%   'ShowViolin'       : Logical. Draw violin plots. Default: true.
%   'YLimMaxSet'       : Override y-axis upper limit.
%   'YLimMinSet'       : Override y-axis lower limit.
%   'Unit'             : String appended to the y-axis label (e.g. 'µm'). Default: ''.
%   'VarNaming'        : Display name for the metric (used when varName is numeric).
%
%  --- ANOVA display ---
%   'PostHocTest'      : Post-hoc method. 'tukey-kramer' (default) | 'tukey' | 'games-howell'.
%   'ANOVAVerbose'     : Print test results to Command Window. Default: true.
%
% OUTPUTS
%   fig         : Figure handle, or [] if no valid data remain after filtering.
%   pVals       : 1-by-nTicks vector of per-condition p-values (currently NaN;
%                 filled by the chosen stat routine via sigPairs).
%   overviewTbl : Table — one row per plotted condition — with columns
%                 Condition, Mean, SD, N (technical reps), NBio (bio-reps).
%   anovaTbl    : Output table from the chosen stat routine, or [] when no
%                 statistical test was run.
%
% DEPENDENCIES
%   Filter_IncludeExclude.m
%   getPlotColors.m
%   getPlotMarkers.m
%   runLMMStat.m  |  runNestedANOVAStat.m  |  runOneWayANOVA.m  |
%   runRepeatedMeasuresANOVA.m  |  runContrastTTests.m

% =========================================================================
% 1. Arguments
% =========================================================================
arguments
    tableT                table
    varName               {mustBeA(varName, ["string","char","double"])}
    options.VarColName    string = ""
    options.VarColIndex   (1,1) double {mustBePositive} = 1
    options.levels        cell = {}

    options.XCat           string = "extrainfoCat"
    options.compareCat     string = "diseaseCat"
    options.compareGroups  cell = {}
    options.ConditionOrder cell = {}

    options.grayband      (1,1) logical = true

    options.includePatterns cell = {}
    options.excludePatterns cell = {}

    options.CutoffBottom  (1,1) double = -Inf
    options.CutoffTop     (1,1) double = Inf
    options.YLimMaxSet    = []
    options.YLimMinSet    = []
    options.YLimMaxForce  (1,1) logical = false   % true = never expand past YLimMaxSet (hide brackets that don't fit)
    options.SplitCat      string = ""
    options.SplitCatLevel string = ""

    options.statTest      string = "welch"
    options.statTest2     string = "lmm"
    options.ANOVA         string = "rm"
    options.ContrastRemove cell = {}
    options.ContrastRemoveG1         cell   = {}
    options.ContrastRemoveG2         cell   = {}
    options.ContrastExcludePatterns  cell   = {}
    options.Contrasts     = []
    options.CollapseRep   (1,1) logical = true
    options.TTestType     string = "paired"

    options.paletteName   string = "distinct"
    options.MarkerSize    (1,1) double {mustBePositive} = 50 % 36
    options.Alpha         (1,1) double {mustBeGreaterThanOrEqual(options.Alpha,0), mustBeLessThanOrEqual(options.Alpha,1)} = 0.75
    options.Jitter        (1,1) double {mustBeNonnegative} = 0.08
    options.ViolinWidth   (1,1) double {mustBePositive} = 0.13
    options.ViolinAlpha   (1,1) double {mustBeGreaterThanOrEqual(options.ViolinAlpha,0), mustBeLessThanOrEqual(options.ViolinAlpha,1)} = 0.25
    options.ShowViolin    (1,1) logical = true

    options.VarNaming     string
    options.Unit          string = ""

    options.MultiplePalette      (1,1) logical = true
    options.CondColors           = []

    options.UseThirdCat          (1,1) logical = true
    options.ThirdCat             string  = "sampletypeCat"
    options.ThirdCatFilter       cell  = {'TS'}

    options.PostHocTest        string = "tukey-kramer"
    options.ANOVAVerbose       (1,1) logical = true

    options.publication        (1,1) logical = true
end

% rename (unchanged block)
levels        = options.levels;
compareCat    = options.compareCat;
xCat          = options.XCat;
compareGroups = string(options.compareGroups);
grayband      = options.grayband;
CutoffBottom  = options.CutoffBottom;
CutoffTop     = options.CutoffTop;
YLimMaxSet    = options.YLimMaxSet;
YLimMinSet    = options.YLimMinSet;
splitCat      = options.SplitCat;
splitCatLevel = options.SplitCatLevel;
% statTest      = options.statTest;
paletteName   = options.paletteName;
mSize         = options.MarkerSize;
alphaVal      = options.Alpha;
jitterAmp     = options.Jitter;
violinWidth   = options.ViolinWidth;
violinAlpha   = options.ViolinAlpha;
showViolin    = options.ShowViolin;
unitName      = options.Unit;

postHocTest          = options.PostHocTest;
anovaVerbose         = options.ANOVAVerbose;
multiplePalette      = options.MultiplePalette;

useThirdCat    = options.UseThirdCat;
thirdCat       = options.ThirdCat;
ThirdCatFilter = options.ThirdCatFilter;
Publication    = options.publication;

if isempty(options.Contrasts)
    options.Contrasts = "all";
end

%outFolder = fileparts(fileparts(tableT.Path{1,1}));

% if varName is not chr
if ischar(varName) || isstring(varName)
    varNameStr = string(varName);
    metricIsName = true;
else
    % numeric vector case — requires VarColName + VarColIndex to re-extract
    varNameStr   = options.VarNaming;
    metricIsName = false;
    if options.VarColName == ""
        error('When passing a numeric vector as varName, you must also supply ''VarColName'' and ''VarColIndex''.');
    end
end

% =========================================================================
% 2. Filter rows via include/exclude on SampleFolder
% =========================================================================
if ~isempty(options.includePatterns) || ~isempty(options.excludePatterns)
    allData = Filter_IncludeExclude(tableT, options.includePatterns, options.excludePatterns);
else
    allData = tableT;
end

% =========================================================================
% 3. Optional pre-filter: keep only rows where splitCat == splitCatLevel
% (old part of function, can be ignored if splitCat is not specified)
% =========================================================================
splitCat      = char(string(splitCat));
splitCatLevel = char(string(splitCatLevel));

if ~isempty(splitCat)
    if ~ismember(splitCat, allData.Properties.VariableNames)
        error('splitCat column "%s" not found in table.', splitCat);
    end
    if isempty(splitCatLevel)
        error('splitCatLevel must be specified when splitCat is used.');
    end
    colSplit = string(allData.(splitCat));
    allData  = allData(colSplit == string(splitCatLevel), :);
    if height(allData) == 0
        warning('No rows remain after splitCat filter "%s" == "%s".', splitCat, splitCatLevel);
        fig = []; pVals = []; overviewTbl = table(); anovaTbl = [];
        return;
    end
end

% Extract metric from the already-filtered allData
if metricIsName
    metricRaw = allData.(varNameStr);
else
    metricRaw = allData.(options.VarColName)(:, options.VarColIndex);
end

% --- Parallel SD column (varName_SD), used for T-marker error bars on mean dots ---
sdColName = varNameStr + "_SD";
if ismember(sdColName, allData.Properties.VariableNames)
    sdRaw = allData.(sdColName);
else
    sdRaw = nan(height(allData), 1);
end

% =========================================================================
% 4. If Levels is empty → infer from table using XCat
% =========================================================================
if isempty(levels)
    if ismember(xCat, allData.Properties.VariableNames)
        levels = reshape(cellstr(string(unique(allData.(xCat)))), 1, []);
    else
        error("Column '%s' not found in table.", xCat);
    end
end

% =========================================================================
% 5. Determine compare categories automatically if not supplied
% =========================================================================
compareVals = unique(string(allData.(compareCat)), 'stable');
compareVals = compareVals(compareVals ~= "");

if isempty(compareGroups)
    compareGroups = compareVals;
end

nGroups = numel(compareGroups);

% =========================================================================
% 6. Apply cutoffs and build filter mask
% =========================================================================
%metric     = double(allData.(varName));
%metricRaw = allData.(varName);

if iscell(metricRaw)
    metric = nan(size(metricRaw)); % preserve size
    for i = 1:numel(metricRaw)
        if ~isempty(metricRaw{i}) && isnumeric(metricRaw{i})
            metric(i) = metricRaw{i};
        else
            metric(i) = NaN;
        end
    end
else
    metric = double(metricRaw);
end

% --- Convert SD column the same way ---
if iscell(sdRaw)
    sdVec = nan(size(sdRaw));
    for i = 1:numel(sdRaw)
        if ~isempty(sdRaw{i}) && isnumeric(sdRaw{i})
            sdVec(i) = sdRaw{i};
        else
            sdVec(i) = NaN;
        end
    end
else
    sdVec = double(sdRaw);
end

metric(metric <  CutoffBottom) = NaN;
metric(metric >  CutoffTop)    = NaN;

compareStr = string(allData.(compareCat));
isEither   = ismember(compareStr, compareGroups);
% isGrp1     = (compareStr == compare1);
% isGrp2     = (compareStr == compare2);
% isEither   = isGrp1 | isGrp2;

xCatCol = string(allData.(xCat));
% xCat filtered to levels
xCatCat = categorical(xCatCol, levels, 'Ordinal', true);
% xCat all for statistics
xCatCatStat = categorical(xCatCol, cellstr(string(unique(xCatCol))), 'Ordinal', false);

%folder  = categorical(string(allData.DataFolder));
folder = categorical(string(allData.expCat));

% valid rows filtered to levels etc.
if useThirdCat
    thirdColAll = string(allData.(thirdCat));
    % filter tp ThirdCatFilter
    ThirdCatCat = categorical(thirdColAll, ThirdCatFilter, 'Ordinal', true);
    % find validrows
    validRows = ~isnan(metric) & isEither & ~ismissing(xCatCat) & ~ismissing(folder) & ~ismissing(ThirdCatCat);
else
    validRows = ~isnan(metric) & isEither & ~ismissing(xCatCat) & ~ismissing(folder);
end
metricF    = metric(validRows);
sdF        = sdVec(validRows);
compareF   = compareStr(validRows);
xCatF      = xCatCat(validRows);
folderF    = folder(validRows);
if isempty(metricF)
    warning('No valid "%s" entries after filtering.', varNameStr);
    fig = []; pVals = []; overviewTbl = table(); anovaTbl = [];
    return;
end

% valid rows for all levels etc. for statistics
validRowsStat = ~isnan(metric) & isEither & ~ismissing(xCatCatStat) & ~ismissing(folder);
metricStat  = metric(validRowsStat);
compareStat = compareStr(validRowsStat);
xCatStat    = xCatCatStat(validRowsStat);
sdStat      = sdVec(validRowsStat);

% --- Third category (optional) ---
if useThirdCat && ismember(thirdCat, allData.Properties.VariableNames)
    %thirdColAll  = string(allData.(thirdCat));
    thirdColF    = thirdColAll(validRows);
    thirdColStat = thirdColAll(validRowsStat);
    thirdLevels  = reshape(cellstr(string(unique(thirdColF))), 1, []);
    % only keep useThirdCat=true if there is actually >1 level or it was explicitly requested
    useThirdCat  = true;
else
    useThirdCat  = false;          % <-- force off; guards the §9 loop branch
    thirdColF    = strings(sum(validRows), 1);
    thirdColStat = strings(sum(validRowsStat), 1);
    thirdLevels  = {};             % <-- empty, not {""}, so §9 loop is never entered
end

% Per-row combo strings (used everywhere instead of prefix-stripping)
if useThirdCat
    comboF    = makeComboStr(string(compareF),    string(xCatF),    thirdColF);
    comboStat = makeComboStr(string(compareStat), string(xCatStat), thirdColStat);
else
    comboF    = makeComboStr(string(compareF),    string(xCatF));
    comboStat = makeComboStr(string(compareStat), string(xCatStat));
end

% =========================================================================
% 7. Levels actually present in filtered data
% =========================================================================
levelsPresent = levels(ismember(levels, cellstr(unique(xCatF))));
K = numel(levelsPresent);
if K == 0
    warning('No xCat levels present after filtering.');
    fig = []; pVals = []; overviewTbl = table(); anovaTbl = [];
    return;
end

% =========================================================================
% 8. Colors and markers
% =========================================================================
folderCats = unique(folderF, 'stable');
nF         = numel(folderCats);
colorsF    = getPlotColors(paletteName, nF);
markersF   = getPlotMarkers(nF);

% Resolve condition color palettes (mirrors Plot_CompareVar logic)
% Build a flat map: 'CT_rest' -> palette name
condColorMap = resolveCondColorMap(compareGroups, levelsPresent, options.CondColors);

% =========================================================================
% 9. Layout constants
% =========================================================================
% Build ordered combined condition list
allCombos = {};
for g = 1:nGroups
    for k = 1:K
        if useThirdCat
            for q = 1:numel(thirdLevels)
                allCombos{end+1} = char(makeComboStr(compareGroups(g), levelsPresent{k}, thirdLevels{q}));
            end
        else
            allCombos{end+1} = char(makeComboStr(compareGroups(g), levelsPresent{k}));
        end
    end
end

% Keep only combinations that actually have data
presentCombos = unique(cellstr(comboStat));   % real combos from data
allCombos = allCombos(ismember(allCombos, presentCombos));

if ~isempty(options.ConditionOrder)
    conditionOrder = cellstr(options.ConditionOrder);
    % keep only those actually present in data
    conditionOrder = conditionOrder(ismember(conditionOrder, allCombos));
    % append any present combos not explicitly listed
    missing = allCombos(~ismember(allCombos, conditionOrder));
    conditionOrder = [conditionOrder, missing];
else
    conditionOrder = allCombos;
end

nTicks   = numel(conditionOrder);
meanHalf = max(jitterAmp, 0.08);
xScale = 0.4;    % decrease to bring ticks closer

% =========================================================================
% 10. Pre-compute per-condition n (for labels)
% =========================================================================
nGrpFlat = zeros(1, nTicks);

for t = 1:nTicks
    nGrpFlat(t) = sum((comboF == conditionOrder{t}) & ~isnan(metricF));
end

pVals = nan(1, nTicks);

% =========================================================================
% 11. Statistics via LMM or nested ANOVA
% =========================================================================
anovaTbl = [];
sigPairs = {};

combinedGroup = cellstr(comboStat);
groupVec      = cellstr(combinedGroup);
uniqueGroups  = unique(groupVec);
nUniqueGroups = numel(uniqueGroups);

% --- build ordered combined levels (preserve your existing x-axis order) ---
% levelsPresent already lists the xCat levels in plot order;
% Build ordered combined levels for stat functions
groupOrder         = unique(cellstr(string(compareStat)), 'stable');
allXCatForStat     = cellstr(string(unique(xCatStat)));
combinedLevelsStat = {};
for g = 1:numel(groupOrder)
    if useThirdCat
        allThirdForStat = cellstr(string(unique(thirdColStat)));
        for k = 1:numel(allXCatForStat)
            for q = 1:numel(allThirdForStat)
                combinedLevelsStat{end+1} = char(makeComboStr(groupOrder{g}, allXCatForStat{k}, allThirdForStat{q}));
            end
        end
    else
        for k = 1:numel(allXCatForStat)
            combinedLevelsStat{end+1} = char(makeComboStr(groupOrder{g}, allXCatForStat{k}));
        end
    end
end
combinedLevelsStat = combinedLevelsStat(ismember(combinedLevelsStat, combinedGroup));
combinedLevelsStat = cellstr(combinedLevelsStat);

% Build a sub-table aligned with the already-filtered data
filteredIdx  = find(validRowsStat);
statSubTbl   = allData(filteredIdx, :);

if strcmp(options.statTest2, 'lmm')
    [sigPairs, anovaTbl, lmeModel, workTbl] = runLMMStat( ...
        statSubTbl, metricStat, ...
        'CombinedVec',     groupVec, ...          
        'BioRepCol',       'expCat', ...
        'ConditionLevels', combinedLevelsStat, ...
        'Contrasts',       options.Contrasts,... % { {'CT_rest','CCM_rest'} }, ... % all
        'ContrastRemove',  options.ContrastRemove,...
        'ContrastExcludePatterns', options.ContrastExcludePatterns,...
        'CollapseRep',      options.CollapseRep,...
        'Correction',      'fdr', ... % fdr
        'Verbose',         anovaVerbose);
        % explicit cell Contrast: { {'CT_rest','CCM_rest'}, {'CT_1kPa','CT_rest'} }
elseif strcmp(options.statTest2, 'nested')
    [sigPairs, anovaTbl, lmeModel, workTbl] = runNestedANOVAStat( ...
        statSubTbl, metricStat, ...
        'ConditionCol',    xCat, ...
        'GroupCol',        compareCat, ...
        'BioRepCol',       'expCat', ...
        'ConditionLevels', levelsPresent, ...
        'GroupLevels',     cellstr(compareGroups), ...
        'Contrasts',       'all', ...       % or 'within' or explicit cell
        'Correction',      'fdr', ...
        'Verbose',         anovaVerbose);
elseif strcmp(options.statTest2, 'onewayTech')
    [sigPairs, anovaTbl, lmeModel, workTbl] = runOneWayANOVA(statSubTbl, metricStat, ...
        'CombinedVec', groupVec, ...
        'MeanBioRep', false, ...
        'Correction', 'fdr');
elseif strcmp(options.statTest2, 'onewayBio')
    [sigPairs, anovaTbl, lmeModel, workTbl] = runOneWayANOVA(statSubTbl, metricStat, ...
        'CombinedVec', groupVec, ...
        'MeanBioRep', true, ...
        'Correction', 'fdr');
elseif strcmp(options.statTest2, 'rm')
    [sigPairs, anovaTbl, lmeModel, workTbl] = runRepeatedMeasuresANOVA(statSubTbl, metricStat, ...
        'CombinedVec', groupVec, ...
        'ConditionLevels', combinedLevelsStat, ...
        'Correction', 'hsd'); % bonferroni
elseif strcmp(options.statTest2, 'ttest')
    [sigPairs, statsTbl, workTbl] = runContrastTTests(statSubTbl, metricStat, ...
        'CombinedVec', groupVec, ...
        'Contrasts',  options.Contrasts,...
        'MeanBioRep', options.CollapseRep, ...
        'TestType',   options.TTestType);
elseif strcmp(options.statTest2, 'wlm')   % weighted linear model on precomputed diffs
    [sigPairs, statsTbl, mdl, workTbl] = runWeightedContrastTest(metricStat, sdStat, groupVec, ...
        'Contrasts',  options.Contrasts, ...
        'Correction', 'fdr', ...
        'Verbose',    anovaVerbose);
    anovaTbl = statsTbl;   % keep output-variable naming consistent with rest of function
end

% =========================================================================
% 12. Draw figure
% =========================================================================
fig = figure('Color', 'w');
ax  = axes('Parent', fig);
hold(ax, 'on');

% =========================================================================
% 12a. Per-combination scatter colors (MultiplePalette)
% =========================================================================
% Order: compare1_lvl1, compare2_lvl1, compare1_lvl2, compare2_lvl2, ...
% scatterColorMap = containers.Map('KeyType','char','ValueType','any');
% 
% if multiplePalette
%     if isempty(multiplePaletteNames)
%         error('MultiplePaletteNames must be provided when MultiplePalette is true.');
%     end
%     comboIdx = 0;
%     for k = 1:K
%         lvl = levelsPresent{k};
%         for g = 1:nGroups
%             cmp = compareGroups(g);
%             comboIdx = comboIdx + 1;
%             key = char(cmp + "_" + string(lvl));
%             if comboIdx <= numel(multiplePaletteNames)
%                 palName = multiplePaletteNames{comboIdx};
%             else
%                 warning('Not enough palette names supplied; recycling last one.');
%                 palName = multiplePaletteNames{end};
%             end
%             scatterColorMap(key) = getPlotColors(palName, nF);  % nF shades
%         end
%     end
% end

% =========================================================================
% 12b. Violin plots
% =========================================================================
if showViolin
    violin1Color = [0.70 0.70 0.70];
    nKde = 256;
    for t = 1:nTicks
        xPos = t * xScale;
        dataG = metricF(comboF == conditionOrder{t});
        dataG = dataG(~isnan(dataG));
        if numel(dataG) >= 3
            drawViolin(ax, dataG, xPos, violinWidth, nKde, ...
                violin1Color, violinAlpha, sprintf('violin_%s', conditionOrder{t}));
        end
    end
end

% Scatter + mean lines
lgHandles = gobjects(nF, 1);
hasHandle = false(nF, 1);

[grpIdx, ~]   = findgroups(comboF, folderF); 
techRepCounts = splitapply(@(v) sum(~isnan(v)), metricF, grpIdx);
showMeanMarkers = any(techRepCounts > 1);

for t = 1:nTicks
    xPos = t * xScale;
    idxG = (comboF == conditionOrder{t});
    if ~any(idxG), continue; end

    yAll   = metricF(idxG);
    xMeans = nan(nF,1);
    yMeans = nan(nF,1);
    sdMeans = nan(nF,1);
    colorsMean = nan(nF,3);
    
    % Extract base condition key, e.g. "CT_rest" from "CT_rest_TS"
    condStr = string(conditionOrder{t});
    parts   = split(condStr, "_");
    baseKey = char(join(parts(1:min(2,end)), "_"));

    for f = 1:nF
        idxF = idxG & (folderF == folderCats(f));
        yy = metricF(idxF);
        if isempty(yy), continue; end

        xx = xPos + jitterAmp * (rand(numel(yy),1) - 0.5);

        % Color: use CondColors palette if provided, otherwise fall back
        % to colorsF which is already built from paletteName
        if multiplePalette && ~isempty(options.CondColors) && isKey(condColorMap, baseKey)
            dotColor = getPlotColors(condColorMap(baseKey), nF);
            dotColor = dotColor(f,:);
        else
            dotColor = colorsF(f,:);
        end

        % if multiplePalette
            % % --- Extract base key (e.g., "CT_rest" from "CT_rest_TL")
            % condStr = string(conditionOrder{t});
            % parts = split(condStr, "_");
            % if numel(parts) >= 2
            %     comboKey = parts(1) + "_" + parts(2);   % -> CT_rest
            % else
            %     comboKey = condStr;
            % end
            % % --- Find matching key in map (partial match)
            % keyList = keys(scatterColorMap);
            % matchIdx = contains(keyList, comboKey);
            % if any(matchIdx)
            %     matchedKey = keyList{find(matchIdx,1,'first')};
            %     comboColors = scatterColorMap(matchedKey);
            %     dotColor = comboColors(f,:);
            % else
            %     dotColor = colorsF(f,:);   % fallback
            % end
            % --- Derive per-folder shade from the condition's palette ---

        % Technical replicate scatter
        h = scatter(ax, xx, yy, mSize, ...
            'Marker',          markersF{f}, ...
            'MarkerFaceColor', dotColor, ...
            'MarkerEdgeColor', dotColor, ...
            'MarkerFaceAlpha', alphaVal, ...
            'MarkerEdgeAlpha', alphaVal, ...
            'DisplayName',     char(folderCats(f)));
        if ~hasHandle(f), lgHandles(f) = h; hasHandle(f) = true; end

        % Mean marker per biological replicate (folder), plotted last
        mF = mean(yy, 'omitnan');
        if isfinite(mF)
            meanOffsets = linspace(-0.04, 0.04, nF);
            xMeans(f) = xPos + meanOffsets(f) + 0.02*(rand-0.5);
            yMeans(f) = mF;
            colorsMean(f,:) = dotColor;
            sdMeans(f) = mean(sdF(idxF), 'omitnan');   % SD is constant per bio-rep/condition; mean() is just a safe reduction
        end
    end

    % Plot mean markers (after technical loop so no overlap)
    for f = 1:nF
        if ~isnan(yMeans(f))
            scatter(ax, xMeans(f), yMeans(f), mSize * 2, ...
                'Marker', markersF{f}, ...
                'MarkerFaceColor', colorsMean(f,:), ...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceAlpha', 1.0, ...
                'MarkerEdgeAlpha', 1.0, ...
                'LineWidth', 1.8, ...
                'HandleVisibility', 'off');
            if ~showMeanMarkers
                % Only 1 tech rep per bio-rep everywhere: mean dot would be
                % redundant with the raw scatter point, so show SD instead
                if isfinite(sdMeans(f)) && sdMeans(f) > 0
                    drawErrorBarT(ax, xMeans(f), yMeans(f), sdMeans(f), 0.018, [0 0 0], 1.2);
                end
            end
        end
    end

    m = mean(yAll, 'omitnan');
    plot(ax, [xPos-meanHalf, xPos+meanHalf], [m m], 'k-', 'LineWidth', 1.8);
end

% =========================================================================
% 13. Fix x-limits and stretch gray bands to full y-extent
% =========================================================================
xlim(ax, [0.5*xScale, (nTicks+0.5)*xScale]);
set(fig, 'Renderer', 'opengl');

validMetric = metricF(~isnan(metricF));
padAll = 0.08 * range(validMetric);
if ~isfinite(padAll) || padAll == 0
    padAll = 0.05 * max(abs(validMetric));
    if ~isfinite(padAll), padAll = 1; end
end
yl = ylim(ax);
if ~isempty(YLimMinSet), yl(1) = YLimMinSet; end
if ~isempty(YLimMaxSet)
    ylim(ax, [yl(1), YLimMaxSet]);      % exact, no padAll
else
    ylim(ax, [yl(1), yl(2) + padAll]);
end
yl = ylim(ax);

if grayband
    bands = findobj(ax, 'Tag', 'band_*');
    for k = 1:K
        b = findobj(ax, 'Tag', sprintf('band_%d', k));
        if isempty(b), continue; end
        set(b, 'YData', [yl(1) yl(1) yl(2) yl(2)]);
    end
    
    uistack(bands, 'bottom');
else
end

if showViolin
    violinPatches = findobj(ax, '-regexp', 'Tag', '^violin');
    if ~isempty(violinPatches), uistack(violinPatches, 'bottom'); end
end

% =========================================================================
% 14 & 15. Significance brackets (one-way ANOVA post-hoc)
% =========================================================================
if ~isempty(sigPairs)
    groupXPos = containers.Map('KeyType','char','ValueType','double');
    for t = 1:nTicks
        groupXPos(conditionOrder{t}) = t * xScale;
    end
    sigPairs = sigPairs(cellfun(@(p) isKey(groupXPos, p{1}) && isKey(groupXPos, p{2}), sigPairs)); % only keep sigPairs where both keys exist

    pairInfo = [];
    % --- in section 14/15, the sigPairs loop ---
    for r = 1:numel(sigPairs)
        gA = sigPairs{r}{1};  gB = sigPairs{r}{2};
        if ~isKey(groupXPos, gA) || ~isKey(groupXPos, gB)
            fprintf('WARNING: skipped pair "%s" vs "%s" — key not in groupXPos\n', gA, gB);
            fprintf('  Available keys: %s\n', strjoin(keys(groupXPos), ', '));
            continue   % <-- was missing; execution fell through to the crash
        end
        xA = groupXPos(gA);   xB = groupXPos(gB);

        pairInfo(end+1).xA   = min(xA, xB);
        pairInfo(end).xB     = max(xA, xB);
        pairInfo(end).star   = sigPairs{r}{3};
        pairInfo(end).span   = abs(xB - xA);
    end

    if ~isempty(pairInfo)
        [~, sortOrd] = sort([pairInfo.span]);
        pairInfo = pairInfo(sortOrd);

        % Initialise per-column ceiling using actual x positions
        xPosList = zeros(1, nTicks);
        colTop   = zeros(1, nTicks);
        
        for i = 1:nTicks
            xPosList(i) = groupXPos(conditionOrder{i});
            vals = metricStat(strcmp(combinedGroup, conditionOrder{i}));
            vals = vals(~isnan(vals));
            if isempty(vals) || ~any(isfinite(vals))
                yl = ylim(ax);
                v  = yl(1);
            else
                v = max(vals(isfinite(vals)));
            end
            colTop(i) = v + 0.04 * (yl(2) - yl(1));
        end

                yl = ylim(ax);
        bracketStep  = 0.07 * (yl(2) - yl(1));
        starPad      = 0.5 * bracketStep;     % headroom for the star above the bar
        capActive    = ~isempty(YLimMaxSet);
        yNeeded      = -Inf;                  % highest point any bracket + star needs
        nSkipped     = 0;
        bracketColor = [0 0 0];

        if options.publication
            bracketLW  = 1.5;
            starFontSize = 16;
        else
            bracketLW  = 1.2;
            starFontSize = 10;
        end

        for r = 1:numel(pairInfo)
            xA   = pairInfo(r).xA;
            xB   = pairInfo(r).xB;
            star = pairInfo(r).star;

            inSpan = xPosList >= xA & xPosList <= xB;
            yBase  = max(colTop(inSpan)) + bracketStep;

            yNeeded = max(yNeeded, yBase + starPad);
            yl = ylim(ax);

            if capActive
                % Fixed upper limit: never touch ylim inside the loop
                if options.YLimMaxForce && (yBase + starPad > YLimMaxSet)
                    nSkipped = nSkipped + 1;
                    continue                    % bracket doesn't fit -> skip it
                end
            else
                % Original behaviour: grow the axis as needed
                if yBase > yl(2)
                    ylim(ax, [yl(1), yBase + bracketStep]);
                    yl = ylim(ax);
                end
            end

            idxA = find(abs(xPosList - xA) < 1e-9, 1);
            idxB = find(abs(xPosList - xB) < 1e-9, 1);

            footGap = 0.025 * (yl(2) - yl(1));
            yFootA  = colTop(idxA) + footGap;
            yFootB  = colTop(idxB) + footGap;

            plot(ax, [xA xA], [yFootA, yBase], 'Color', bracketColor, 'LineWidth', bracketLW);
            plot(ax, [xB xB], [yFootB, yBase], 'Color', bracketColor, 'LineWidth', bracketLW);
            plot(ax, [xA xB], [yBase, yBase],  'Color', bracketColor, 'LineWidth', bracketLW);
            text(ax, (xA+xB)/2, yBase + 0.005*(yl(2)-yl(1)), star, ...
                'HorizontalAlignment', 'center', 'FontSize', starFontSize, 'Color', bracketColor);

            colTop(inSpan) = yBase;
        end

        % --- final y-limit decision ---
        yl = ylim(ax);
        if capActive
            if options.YLimMaxForce
                if nSkipped > 0
                    warning('%s: %d significant bracket(s) hidden. YLimMaxSet >= %.3g is needed (currently %.3g).', ...
                        varNameStr, nSkipped, yNeeded, YLimMaxSet);
                end
            elseif yNeeded > YLimMaxSet
                warning('%s: YLimMaxSet = %.3g is too low for the brackets; extended to %.3g.', ...
                    varNameStr, YLimMaxSet, yNeeded);
                ylim(ax, [yl(1), yNeeded]);
            end
        else
            ylim(ax, [yl(1), yl(2) + 0.04*(yl(2)-yl(1))]);
        end
        % yl = ylim(ax);
        % bracketStep  = 0.07 * (yl(2) - yl(1));
        % bracketColor = [0 0 0];
        % 
        % % Choose line width and font size based on publication mode
        % if options.publication
        %     bracketLW  = 1.5;
        %     starFontSize = 16;
        % else
        %     bracketLW  = 1.2;
        %     starFontSize = 10;
        % end
        % 
        % for r = 1:numel(pairInfo)
        %     xA   = pairInfo(r).xA;
        %     xB   = pairInfo(r).xB;
        %     star = pairInfo(r).star;
        % 
        %     inSpan = xPosList >= xA & xPosList <= xB;
        %     yBase  = max(colTop(inSpan)) + bracketStep;
        % 
        %     yl = ylim(ax);
        %     if yBase > yl(2)
        %         ylim(ax, [yl(1), yBase + bracketStep]);
        %         yl = ylim(ax);
        %     end
        % 
        %     % Foot of each vertical line = current colTop for that specific column
        %     % (i.e. the top of whatever bracket or data was last drawn there)
        %     idxA = find(abs(xPosList - xA) < 1e-9, 1);
        %     idxB = find(abs(xPosList - xB) < 1e-9, 1);
        % 
        %     footGap = 0.025 * (yl(2) - yl(1));
        %     yFootA  = colTop(idxA) + footGap;
        %     yFootB  = colTop(idxB) + footGap;
        % 
        %     % vertical lines from previous bracket top (staggered) up to new horizontal bar
        %     plot(ax, [xA xA], [yFootA, yBase], 'Color', bracketColor, 'LineWidth', bracketLW);
        %     plot(ax, [xB xB], [yFootB, yBase], 'Color', bracketColor, 'LineWidth', bracketLW);
        %     % horizontal bar
        %     plot(ax, [xA xB], [yBase, yBase],  'Color', bracketColor, 'LineWidth', bracketLW);
        %     % significance label
        %     text(ax, (xA+xB)/2, yBase + 0.005*(yl(2)-yl(1)), star, ...
        %         'HorizontalAlignment', 'center', 'FontSize', starFontSize, 'Color', bracketColor);
        % 
        %     % update colTop for all columns spanned by this bracket
        %     colTop(inSpan) = yBase;
        % end
        % 
        % yl = ylim(ax);
        % ylim(ax, [yl(1), yl(2) + 0.04*(yl(2)-yl(1))]);
    end
end

% =========================================================================
% 16. N labels
% =========================================================================
if ~options.publication
    yl = ylim(ax);
    nLabelPad = 0.01 * (yl(2) - yl(1));
    
    for t = 1:nTicks
        xPos = t * xScale;
        yl   = ylim(ax);
        yTop = yl(2) + nLabelPad;
        if nGrpFlat(t) > 0
            text(ax, xPos, yTop, sprintf('n=%d', nGrpFlat(t)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment',   'bottom', ...
                'FontSize', 7, 'Color', [0.3 0.3 0.3], 'Clipping', 'off');
        end
    end
    
    if isempty(YLimMaxSet)
        yl = ylim(ax);
        ylim(ax, [yl(1), yl(2) + 0.05 * (yl(2) - yl(1))]);
    end
end

% =========================================================================
% 17. Axes decoration
% =========================================================================
xticks(ax, (1:nTicks) * xScale);
xticklabels(ax, conditionOrder);
ax.TickLabelInterpreter = 'none';
ylabel(ax, sprintf('%s %s', varNameStr, unitName), 'Interpreter', 'none');

if ~options.publication
    titleStr = sprintf('%s  |  %s', varNameStr, strjoin(compareGroups, " vs "));
    title(ax, titleStr, 'Interpreter', 'none');
    
    testLabel = sprintf('Stat = %s', options.statTest2);
    
    if grayband
        subtitleParts = {sprintf('%s | gray band = %s', testLabel, compareGroups(end))};
    else
        subtitleParts = {sprintf('%s', testLabel)};
    end
    
    % if grayband
    %     subtitleParts = {sprintf('one-way ANOVA, post-hoc: %s | gray band = %s', postHocTest, compareGroups(end))};
    % else
    %     subtitleParts = {sprintf('one-way ANOVA, post-hoc: %s', postHocTest)};
    % end
    if isfinite(CutoffBottom), subtitleParts{end+1} = sprintf('cut_bot=%.4g', CutoffBottom); end
    if isfinite(CutoffTop),    subtitleParts{end+1} = sprintf('cut_top=%.4g', CutoffTop);    end
    subtitle(ax, strjoin(subtitleParts, '  |  '), 'FontSize', 8, 'Color', [0.4 0.4 0.4]);
end

grid(ax, 'on');
ax.Layer = 'bottom';  % push grid lines behind all plotted content
set(ax, 'Box', 'off');
ax.Layer = 'bottom';
% ax.XGrid = 'off'; % turn vertical grid line off
% ax.YGrid = 'on';
grid(ax, 'minor');
ax.XMinorGrid = 'off';
%ax.YMinorGrid = 'on';

% --- Publication styling ---
if options.publication
    pubFontSize  = 14;
    pubLineWidth = 1.5;
    set(ax, ...
        'FontSize',        pubFontSize, ...   % all text in axis
        'LineWidth',       pubLineWidth, ...  % thickness of box around graph
        'TickLength',      [0.015 0.025], ... % thick mark length
        'TickDir',         'in',...           % direction thicks (in/out)
        'GridLineWidth',   pubLineWidth*0.5, ...  % major grid line thickness
        'GridAlpha',       0.1, ...           % major grid opacity
        'MinorGridLineWidth', pubLineWidth * 0.25, ...  % minor grid line thickness
        'MinorGridAlpha',  0.2);              % minor grid opacity
    ax.YLabel.FontSize = pubFontSize*1.1;     % fontsize for ylabel, xlabel
    ax.XAxis.LineWidth = pubLineWidth;        % x-asis line thickness
    ax.YAxis.LineWidth = pubLineWidth;        % y-asis line thickness
end

useIdx = find(hasHandle);
labels = cellstr(folderCats(useIdx));
labelsShort = cellfun(@(s) shortenLabel3(s), labels, 'UniformOutput', false);

function out = shortenLabel3(s)
    parts = strings(1,3);
    rest = char(s);
    for i = 1:3
        [tok, rest] = strtok(rest, '_');
        if isempty(tok)
            parts = parts(1:i-1);
            break
        end
        parts(i) = string(tok);

        if isempty(rest)
            parts = parts(1:i);
            break
        end
        rest = rest(2:end); % remove leading underscore
    end
    out = strjoin(parts, '_');
end

if ~isempty(useIdx) && ~options.publication
    legend(lgHandles(useIdx), labelsShort, ...
        'Location', 'southoutside', 'Interpreter', 'none', ...
        'NumColumns', min(nF, 3));
end

% =========================================================================
% 18. Overview table
% =========================================================================
overviewTbl = table( ...
    conditionOrder(:), ...
    nan(nTicks,1), nan(nTicks,1), nan(nTicks,1), nan(nTicks,1), ...
    'VariableNames', {'Condition','Mean','SD','N','NBio'});

for t = 1:nTicks
    idx = (comboF == conditionOrder{t});
    data = metricF(idx);
    data = data(~isnan(data));

    overviewTbl.Mean(t) = mean(data, 'omitnan');
    overviewTbl.SD(t)   = std(data,  'omitnan');
    overviewTbl.N(t)    = numel(data);

    filteredIdx = find(validRows);
    subTbl = allData(filteredIdx(idx), :);
    if ismember('expCat', subTbl.Properties.VariableNames)
        overviewTbl.NBio(t) = numel(unique(subTbl.expCat));
    else
        overviewTbl.NBio(t) = NaN;
    end
end

fprintf('\n--- Overview: %s ---\n', strjoin(conditionOrder, ' | '));
disp(overviewTbl);

hold(ax, 'off');

% SAVING
% exportgraphics(fig, fullfile(outFolder, sprintf('%s_plot.png', varName)), 'Resolution', 300);
% if runANOVA
%     writetable(overviewTbl, fullfile(outFolder, sprintf('%s_overview.xlsx', varName)));
%     if iscell(anovaTbl)
%         anovaTbl = cell2table(anovaTbl(2:end,:), ...
%             'VariableNames', string(anovaTbl(1,:)));
%     end
%     writetable(anovaTbl, fullfile(outFolder, sprintf('%s_anova.xlsx', varName)));
% end
end


% =========================================================================
% LOCAL HELPER: drawViolin (unchanged)
% =========================================================================
% function out = pToStars(p)
%     if p < 0.001,     out = '***';
%     elseif p < 0.01,  out = '**';
%     elseif p < 0.05,  out = '*';
%     else,             out = 'ns';
%     end
% end

function drawViolin(ax, data, xCenter, halfWidth, nKde, faceColor, faceAlpha, tagStr)
    if numel(data) < 3, return; end

    dataMin = min(data);  dataMax = max(data);
    span    = dataMax - dataMin;
    if span == 0
        yq  = dataMin + [-1e-9; 1e-9];
        den = [1; 1];
    else
        pad = 0.10 * span;
        yq  = linspace(dataMin - pad, dataMax + pad, nKde)';
        den = ksdensity(data, yq);
    end

    maxDen = max(den);
    if maxDen == 0, return; end
    denNorm = (den / maxDen) * halfWidth;

    xRight = xCenter + denNorm;
    xLeft  = xCenter - denNorm;
    xPoly  = [xRight; flipud(xLeft)];
    yPoly  = [yq;     flipud(yq)];

    patch(ax, xPoly, yPoly, faceColor, ...
        'FaceAlpha', faceAlpha, ...
        'EdgeColor', faceColor * 0.6, ...
        'EdgeAlpha', min(faceAlpha + 0.3, 1), ...
        'LineWidth', 0.8, ...
        'Tag',       tagStr);
end

% -- helper --
function comboStr = makeComboStr(varargin)
% Concatenate 2 or 3 string vectors with "_", skipping empty components.
    result = string(varargin{1});
    for k = 2:numel(varargin)
        part = string(varargin{k});
        % only append if the part is non-empty for all elements
        if any(part ~= "")
            result = result + "_" + part;
        end
    end
    comboStr = result;
end

% -- helper --
function cmap = resolveCondColorMap(compareGroups, levelsPresent, userColors)
% Build containers.Map: 'CT_rest' -> palette string
% If userColors is provided, match keys by substring (case-insensitive),
% not exact match. First match wins.

defaultPalettes = {'blues','oranges','greens','purples', ...
                   'indigos','magentas','teals','redoranges'};
cmap       = containers.Map('KeyType','char','ValueType','char');
paletteIdx = 0;

for g = 1:numel(compareGroups)
    cg = char(compareGroups(g));
    for k = 1:numel(levelsPresent)
        lv  = char(levelsPresent{k});
        key = [cg, '_', lv];

        matchedPalette = '';

        if ~isempty(userColors) && isa(userColors, 'containers.Map')
            userKeys = keys(userColors);
            for u = 1:numel(userKeys)
                if contains(key, userKeys{u}, 'IgnoreCase', true)
                    matchedPalette = userColors(userKeys{u});
                    break   % first match wins
                end
            end
        end

        if ~isempty(matchedPalette)
            cmap(key) = matchedPalette;
        else
            paletteIdx = paletteIdx + 1;
            cmap(key)  = defaultPalettes{mod(paletteIdx-1, ...
                             numel(defaultPalettes)) + 1};
        end
    end
end
end

% helper
function drawErrorBarT(ax, xCenter, yCenter, sd, capHalfWidth, color, lineWidth)
% Draws a vertical T-shaped error bar (mean ± SD) centered at (xCenter, yCenter).
    yTop = yCenter + sd;
    yBot = yCenter - sd;

    % vertical stem
    plot(ax, [xCenter xCenter], [yBot yTop], 'Color', color, ...
        'LineWidth', lineWidth, 'HandleVisibility', 'off');
    % top cap
    plot(ax, [xCenter-capHalfWidth, xCenter+capHalfWidth], [yTop yTop], ...
        'Color', color, 'LineWidth', lineWidth, 'HandleVisibility', 'off');
    % bottom cap
    plot(ax, [xCenter-capHalfWidth, xCenter+capHalfWidth], [yBot yBot], ...
        'Color', color, 'LineWidth', lineWidth, 'HandleVisibility', 'off');
end