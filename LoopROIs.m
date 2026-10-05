%% LOOP analysis for ROIs

% To collapse all the sections press these three keys; Ctrl + ,
close all; clear all;

% User INPUTS to get the data
% Black BG, minSize, folder names !!!

% Pathinfo
MainFolder = 'D:\3_FRETTS_HBMVEC_analysis'; % main directoryq
DateFolders = {'2025_09_22_TLsAndVEFL_ok'}; %{'2025_08_18_StableVECadVin_fix72h_ok'}; % subfolders per experimentday
SampleFolders = {}; % If empty {} -> runs through all the folders in DateFolders

% Optional filtering keywords
includeKeywords = {'VECadTL_stHP35_CT'};    % -> keep only these (needs to contain ALL Keywords, otherwise you switch & to | in line 76)
excludeKeywords = {};                    % -> remove these

% Naming
minSize_Vin = 8;
minSize_VE = 30;
%minSize = 8; % 25-40 for VECad % 8 (no, 15) for Vin
outputPath = '';
dataType = '.tif';
sampleName = 'Cell';
separator = '_';
segmentfolder = 'SegmentD_4'; %SegmentD_4
iTag = 'ch0';
aTag = '_A'; % Tag for acceptor image
data.delIdx = []; %enter the number of the FA that you want to exclude, such as [3 6 10]
cellIdx = [];

flags.ROIDorA               = 1; % 0 = masking on acceptor | 1 = masking on donor (ROI and FAMask process)
flags.manualROI             = 1; % 0 = off | 1 = removal | 2 = keep these
flags.additionalmanualROI   = 1; % 0 = normal behaviour | 1 = load existing ROI and add to it

flags.flagOneFRETimage      = 1; % shows only unadjusted FRET image
flags.AutoFASegments        = 2; % 0 = slider | 1 = Automatic | 2 = tiff from Ilastik (check D or A!)
flags.AutoCellSegment       = 1; % 1 = code | 2 = tiff from Ilastik 
% For acceptor analysis (assumes same datatype & sampleName as donor)
flags.AcceptorData.switch    = 1; % use acceptor data (0: only load donor)

% Load data
% Loop through each date folder
for D = 1:numel(DateFolders)
    currentDatePath = fullfile(MainFolder, DateFolders{D});
    ConditionFolders = SampleFolders;
    if isempty(SampleFolders)
        % Get all subfolders in this date folder (excluding '.' and '..')
        subInfo = dir(currentDatePath);
        ConditionFolders = {subInfo([subInfo.isdir]).name};
        ConditionFolders = ConditionFolders(~ismember(ConditionFolders, {'.', '..'}));
    else
    end

    % ---- FILTERING STEP ----
    % --- Inclusion filtering (if specified)
    if ~isempty(includeKeywords)
        includeMask = true(1, numel(ConditionFolders));
        for k = 1:numel(includeKeywords)
            includeMask = includeMask & contains(ConditionFolders, includeKeywords{k}, 'IgnoreCase', true);
        end
        ConditionFolders = ConditionFolders(includeMask);
    end

    % --- Exclusion filtering (if specified)
    if ~isempty(excludeKeywords)
        excludeMask = false(size(ConditionFolders));
        for k = 1:numel(excludeKeywords)
            excludeMask = excludeMask | contains(ConditionFolders, excludeKeywords{k}, 'IgnoreCase', true);
        end
        ConditionFolders = ConditionFolders(~excludeMask);
    end

    % Loop through each condition folder
    for E = 1:numel(ConditionFolders)
        dataPath = fullfile(currentDatePath, ConditionFolders{E});
        maskPath = dataPath;
        if contains(ConditionFolders{E}, 'Vin')
            minSize = minSize_Vin;
        else
            minSize = minSize_VE;
        end
        cellIdx = []; % also necessary for bigloop

        % Running this folder
        fprintf('Processing folder: %s\n', dataPath);
         
        %% Detect available cells if user left cellIdx empty
        files = dir([dataPath filesep '*' dataType]);
        
        if isempty(cellIdx)
            % Find all files that match the naming pattern for channel 0 (intensity)
            intFiles = dir([dataPath filesep sampleName '*' separator iTag dataType]);
        
            % Extract the cell numbers by parsing file names
            cellNumbers = [];
           
            for f = 1:length(intFiles)
                fname = intFiles(f).name;
                sepPos = strfind(fname, separator); % ADJUSTED
                numStr = extractBetween(fname, length(sampleName)+1, sepPos(1)-1); %ADJUSTED
                cellNumbers(end+1) = str2double(numStr);
            end
        
            % Use sorted unique numbers
            cellIdx = unique(cellNumbers);
        end

        % cellIdx = [14];

        countLoop = 1; %DONT CHANGE THIS! 
        for i = cellIdx
            %% Get Donor images
            intIdx = find(contains({files.name}, [sampleName num2str(i) separator iTag]));
            intName = files(intIdx).name;
            data.int = imread([dataPath filesep intName]);

            %% Get Ilastik mask
            % Get list of possible files
            searchStr = [sampleName num2str(i) separator '*' iTag '*_FASegment.tiff'];
            disp('Constructed mask search string:')
            disp(searchStr)
            
            % Find matching mask file
            maskFiles = dir(fullfile(maskPath, segmentfolder, searchStr)); % Find matching mask file

            % Display what we found (optional)
            if isempty(maskFiles)
                fprintf('No mask found for %s\n', searchStr);
                data.useExternalMask = 0;
            else
                maskName = maskFiles(1).name;  % take first match
                fprintf('Using mask: %s\n', maskName);
                data.mask = imread(fullfile(maskPath, segmentfolder, maskName));
                %data.mask = ~data.mask;  % If BG is white !!!!!!!!!
                data.useExternalMask = 1;
            end

            flags.useExternalMask       = data.useExternalMask; % uses Ilastik mask if provided in the FRET path
        
            % get name
            str = [separator iTag];
            idx = strfind(intName, str);
            intName(idx:idx+length(str)-1) = [];
        
            str = dataType;
            idx = strfind(intName, str);
            intName(idx:idx+length(str)-1) = [];    
            data.name = intName;
        
            % parse
            data.dataPath = dataPath;
            data.outputPath = outputPath;
            data.sampleIndex = [sampleName num2str(i)];
            data.index = i;
            data.logicalSize = 700; 
            
            % LOOP - THIS IS WHERE THE MAGIC HAPPENS
            close all   
            %data.mask = ~data.mask; % OFF for VECad % ON for Vin
            [mask] = ROIandSeg(data, flags, minSize);
            countLoop = countLoop + 1;
        end
        close all;
    end
end   
            
disp('ROIs completed')
close all;