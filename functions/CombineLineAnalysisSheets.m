function CombineLineAnalysisSheets(outputPath)
% Combines all FRET and Intensity sheets from a Line_Analysis Excel file
% in the input folder into Collected_FRET and Collected_Intensity sheets.
inputFolder = outputPath;
emptyRows = 5;

% Find the Excel file containing 'Line_Analysis'
files = dir(fullfile(inputFolder, '*.xlsx'));
lineAnalysisFile = '';
for k = 1:length(files)
    if contains(files(k).name, 'Line_Analysis', 'IgnoreCase', true)
        lineAnalysisFile = fullfile(inputFolder, files(k).name);
        break
    end
end
if isempty(lineAnalysisFile)
    disp('No Line_Analysis Excel file found in the folder. Skipping sheet combination.');
    return
end

inputFile = lineAnalysisFile;
outputFile = fullfile(inputFolder, 'collected_output.xlsx');

% Get sheet names
[~, sheets] = xlsfinfo(inputFile);

% Initialize containers and flags
collectedFRET = {};
collectedIntensity = {};
firstFRETSheet = true;
firstIntensitySheet = true;

% Loop over sheets
for i = 1:length(sheets)
    sheet = sheets{i};
    data = readcell(inputFile, 'Sheet', sheet);
    if isempty(data), continue; end
    
    % Find last non-empty row manually
    lastRow = 0;
    for r = 1:size(data,1)
        row = data(r,:);
        rowHasData = false;
        for c = 1:size(row,2)
            cellValue = row{c};
            if isnumeric(cellValue) && ~isnan(cellValue)
                rowHasData = true; break
            elseif ischar(cellValue) && ~all(isspace(cellValue))
                rowHasData = true; break
            elseif isstring(cellValue) && strlength(cellValue) > 0
                rowHasData = true; break
            end
        end
        if rowHasData, lastRow = r; end
    end
    if lastRow == 0, continue; end
    data = data(1:lastRow,:);
    
    % Replace missing with empty string
    for r = 1:size(data,1)
        for c = 1:size(data,2)
            if ismissing(data{r,c}), data{r,c} = ''; end
        end
    end
    
    % Append to correct collection
    if contains(sheet,'FRET','IgnoreCase',true)
        if ~firstFRETSheet
            data = data(2:end,:);
        else
            firstFRETSheet = false;
        end
        if ~isempty(collectedFRET)
            collectedFRET = [collectedFRET; cell(emptyRows, size(data,2))];
        end
        collectedFRET = [collectedFRET; data];
    elseif contains(sheet,'Intensity','IgnoreCase',true)
        if ~firstIntensitySheet
            data = data(2:end,:);
        else
            firstIntensitySheet = false;
        end
        if ~isempty(collectedIntensity)
            collectedIntensity = [collectedIntensity; cell(emptyRows, size(data,2))];
        end
        collectedIntensity = [collectedIntensity; data];
    end
end

% Write to Excel
if ~isempty(collectedFRET)
    writecell(collectedFRET, outputFile, 'Sheet', 'Collected_FRET');
end
if ~isempty(collectedIntensity)
    writecell(collectedIntensity, outputFile, 'Sheet', 'Collected_Intensity');
end

disp('Line Analysis sheets combined successfully!');
end