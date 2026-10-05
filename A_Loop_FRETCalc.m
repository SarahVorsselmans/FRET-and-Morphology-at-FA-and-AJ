% To collapse all the sections press these three keys; Ctrl + ,
close all; clear all;
%%
run('C:\Program Files\DIPimage 2.9\dipstart.m');

%% ================= 0_User Configuration =================
% Pathinfo
MainFolder = 'D:\3_FRETTS_HBMVEC_analysis'; % main directory
DateFolders = {'2026_04_272727_VinVEGel_20260124_ok'}; % subfolders per experimentday %'2026_01_252627_WellRockGel', 
SampleFolders = {'20260427_VinTS_FL_CT_glass_si72'}; % If empty {} -> runs through all the folders in DateFolders

% Optional filtering keywords
includeKeywords = {};  % 'Clover'        % -> keep only these
excludeKeywords = {};          % -> remove these

% Naming
outputPath = '';
dataType = '.tif';
sampleName = 'Cell';
separator = '_';
segmentfolder = 'SegmentD_4'; %SegmentD_4 %CHANGE FOR VE
iTag = 'ch0';
gTag = 'ch1';
sTag = 'ch2';
aTag = '_A'; % Tag for acceptor image
data.delIdx = []; %enter the number of the FA that you want to exclude, such as [3 6 10]
cellIdx = [];

%% ================= 0_Run Settings =================
% to do: put flags in user interface
flags.LineAnalysis          = 0; %for eccentricity, circularity, orientation, and sub-FA FRET analysis %1
flags.flagSaveMasks         = 0; % saves the masks for the cell, the FA, the cytosol and the background %1
flags.ROIDorA               = 1; % 0 = masking on acceptor | 1 = masking on donor (ROI and FAMask process)
flags.manualROI             = 0; % 0 = off | 1 = removal | 2 = keep these
flags.additionalmanualROI   = 0; % 0 = normal behaviour | 1 = load existing ROI and add to it
flags.ContinueWithoutROI    = 0; % 0 = don't continue and ask for ROI | 1 = continue without ROI

flags.flagTesting           = 0; % for the FA analysis, you dont want to see this.
flags.flagPhasorHistogram   = 0; % shows (and saves) the figure with phasor and histogram %1
flags.flagFRETimage         = 1; % shows (and saves) the FRET images with and without adapted range and intensity overlay
flags.flagOneFRETimage      = 1; % shows only unadjusted FRET image
flags.flagShowFaList        = 0; % shows the the figure of FAs with assigned colors
flags.flagShowFAAreaInt     = 0; % shows graphs of Area vs Total Intensity/Average Intensity %1
flags.flagShowFAFretInt     = 0; % shows graphs of FRET Index vs FA Size/FA Intensity %1
flags.flagFAindex           = 0; %to find the index assigned to a selected FA
flags.CytosolBackground     = 0; % do analysis for cytosol and background
flags.CytosolClean          = 0; % remove background from cytosol (remove very low intenisty pixels; check segmentation and adjust threshold value if needed)
flags.FAMoreAnalysis        = 0; % 0 = off, 1 = on

flags.AutoFASegments        = 2; % 0 = slider | 1 = Automatic | 2 = tiff from Ilastik (check D or A!)
flags.CellSegment           = 0; % 0 = off | 1 = on
flags.AutoCellSegment       = 1; % 1 = code | 2 = tiff from Ilastik 
% For acceptor analysis (assumes same datatype & sampleName as donor)
flags.AcceptorData.switch   = 1; % use acceptor data (0: only load donor)

%% ================= 1_Load Data =================
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
        includeMask = false(size(ConditionFolders));
        for k = 1:numel(includeKeywords)
            includeMask = includeMask | contains(ConditionFolders, includeKeywords{k}, 'IgnoreCase', true);
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
        data.ConditionFolder = ConditionFolders{E};
        clear allOutputMetrics % necessary somehow
        cellIdx = []; % also necessary for bigloop
        clear data.CellSegThreshSwitch % to be safe
        clear data.CellSegThres % to be safe

        % Running this folder
        fprintf('Processing folder: %s\n', dataPath);

        % Load existing results if available % COMMENT THIS OUT IF YOU WANT TO OVERWRITE THE WHOLE FILE
        if exist(append(dataPath, filesep, 'FALCON Output', filesep, 'allOutputMetrics.mat'), 'file')
            loaded = load(append(dataPath, filesep, 'FALCON Output', filesep, 'allOutputMetrics.mat')); % load into struct
            allOutputMetrics = loaded.allOutputMetrics;  % Extract the actual variable
        end

        Info.FALCONpath = fullfile(dataPath, 'FALCON Output');
        if ~exist(Info.FALCONpath, 'dir')
            mkdir(Info.FALCONpath);
        end
        
        % For cell segmentation
        if isfile([dataPath filesep 'CellSegThresh' '.txt'])
             data.CellSegThreshSwitch = 1;
             data.CellSegThresh = readmatrix([dataPath filesep 'CellSegThresh' '.txt']);
        else % if the text file doesn't exist
             data.CellSegThreshSwitch = 0;
        end
        
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
            flimGIdx = find(contains({files.name}, [sampleName num2str(i) separator gTag]));
            flimSIdx = find(contains({files.name}, [sampleName num2str(i) separator sTag]));
            
            intName = files(intIdx).name;
            flimGName = files(flimGIdx).name;
            flimSName = files(flimSIdx).name;
            
            data.int = imread([dataPath filesep intName]);
            data.FLIMG = double(imread([dataPath filesep flimGName]));
            data.FLIMS = double(imread([dataPath filesep flimSName]));
        
            %% Open acceptor
            if flags.AcceptorData.switch
                aIdx = find(contains({files.name}, [sampleName num2str(i) aTag])); % separator not used consistently
                aName = files(aIdx).name; % ch0
                data.acceptor = imread([dataPath filesep aName]);
                data.flags.AcceptorData.switch = 1;
            else
                data.flags.AcceptorData.switch = 0;
            end
            
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
                %data.mask = ~data.mask;  % If BG is white (Turn on for Vin)
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
            
            %get the date that data is obtained (not analysed) 
            dateMatch = regexp(dataPath, '\d{8}', 'match');
        
            if ~isempty(dateMatch)
            data.SampleDate = dateMatch{1};   % store as text (char)
            else
            data.SampleDate = '';             % fallback if nothing found
            end
        
            % parse
            data.dataPath = dataPath;
            data.outputPath = fullfile(data.dataPath, 'FALCON Output');
            data.sampleIndex = [sampleName num2str(i)];
            data.index = i;
            data.logicalSize = 700; 
            
            % LOOP - THIS IS WHERE THE MAGIC HAPPENS
            close all   
            if flags.LineAnalysis
            [Int, FRET_Index_Data, mask, outputMetrics, FA] = A_FRETCalc_VinVECad(data, flags);
            else
            [Int, FRET_Index_Data, mask, outputMetrics] = A_FRETCalc_VinVECad(data, flags);
            end
        
            % Extract data in MATLAB file
            allOutputMetrics(i) = outputMetrics;
            save(fullfile(data.outputPath, 'allOutputMetrics.mat'), 'allOutputMetrics');
           
            % Save to excel file %
            % %parse
            % tmp = strsplit(dataPath, '\');
            % pointers.DirName = tmp{end};
            % pointers.FileName = data.name;
            % pointers.directory = data.dataPath;
            % %save
            % Save2Excel(outputMetrics, pointers, countLoop, flags);
            % if flags.LineAnalysis
            % LineAnalysis2Excel(data, FA, pointers);
            % end

            countLoop = countLoop + 1;
        end
            %CombineLineAnalysisSheets (outputPath);
        close all;
    end
end   
            
disp('FALCON ANALYSIS IS SUCCESSFULLY COMPLETED :)')
close all;
