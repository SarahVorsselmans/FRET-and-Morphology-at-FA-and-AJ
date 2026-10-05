close all; clear all;

%% User Inputs

% Pathinfo
MainFolder = 'D:\3_FRETTS_HBMVEC_analysis'; % main directory
DateFolder = {'2026_03_2222_FRETint'}; %{'2026_01_252627_WellRockGel'}; % subfolders per experimentday
SampleFolder = {'20260322_VinTS_FL_CT_si72'}; % If empty {} -> runs through all the folders in DateFolders

% Optional filtering keywords
includeKeywords = {'VinT'};          % -> keep only these (needs to contain all Keywords, otherwise you switch & to | in line 76)
excludeKeywords = {};          % -> remove these

% Naming
user.Intensity_Threshold = 5; %just to remove very dim pixels % was 8
user.FRET_Max = 50;
user.Image_Resolution = 512;
user.Image_ImSize = 123.2577; % µm
user.Image_PixSize = 0.24073767; % µm
user.minSize_Vin = 8;
user.minSize_VE = 25;  
user.medfiltValue = 9; 

outputPath = '';
dataType = '.tif';
sampleName = 'Cell';
separator = '_';
segmentfolder = 'SegmentD_4'; %SegmentD_4
DDTag = 'ch00'; DATag = 'ch01';
ADTag = 'ch02'; AATag = 'ch03';

flags.manualROI             = 1; % 0 = off | 1 = removal | 2 = keep these

%% Prepare output
allOutputMetrics = table();

%% Load Data
% Loop through each date folder to get all Sample(Condition)Folders
for D = 1:numel(DateFolder)
    currentDatePath = fullfile(MainFolder, DateFolder{D});
    if isempty(SampleFolder)
        % Get all subfolders in this date folder (excluding '.' and '..')
        subInfo = dir(currentDatePath);
        SampleFolder = {subInfo([subInfo.isdir]).name};
        SampleFolder = SampleFolder(~ismember(SampleFolder, {'.', '..'}));
    else
    end

    % ---- FILTERING STEP ----
    %user.SampleFolder = Filter_IncludeExclude(user, includePatterns, excludePatterns);
    % --- Inclusion filtering (if specified)
    if ~isempty(includeKeywords)
        includeMask = true(1, numel(SampleFolder));
        for k = 1:numel(includeKeywords)
            includeMask = includeMask & contains(SampleFolder, includeKeywords{k}, 'IgnoreCase', true);
        end
        SampleFolder = SampleFolder(includeMask);
    end

    % --- Exclusion filtering (if specified)
    if ~isempty(excludeKeywords)
        excludeMask = false(size(SampleFolder));
        for k = 1:numel(excludeKeywords)
            excludeMask = excludeMask | contains(SampleFolder, excludeKeywords{k}, 'IgnoreCase', true);
        end
        SampleFolder = SampleFolder(~excludeMask);
    end

    % Loop through each sample(condition) folder
    for E = 1:numel(SampleFolder)
        dataPath = fullfile(currentDatePath, SampleFolder{E});
        maskPath = dataPath;
        cellIdx = []; % also necessary for bigloop
        outputPath = fullfile(dataPath, 'FALCON Output');
        if ~exist(outputPath, 'dir')
            mkdir(outputPath);
        end

        % Running this folder
        fprintf('Processing folder: %s\n', dataPath);
         
        % Detect available cells if user left cellIdx empty
        files = dir([dataPath filesep '*' dataType]);
        
        if isempty(cellIdx)
            % Find all files that match the naming pattern for channel x
            DDFiles = dir([dataPath filesep sampleName '*' separator DDTag dataType]);
            DAFiles = dir([dataPath filesep sampleName '*' separator DATag dataType]);
            ADFiles = dir([dataPath filesep sampleName '*' separator ADTag dataType]);
            AAFiles = dir([dataPath filesep sampleName '*' separator AATag dataType]);
        
            % Extract the cell numbers by parsing file names
            cellNumbers = [];
           
            for f = 1:length(DDFiles)
                fname = DDFiles(f).name;
                sepPos = strfind(fname, separator);
                numStr = extractBetween(fname, length(sampleName)+1, sepPos(1)-1);
                cellNumbers(end+1) = str2double(numStr);
            end
        
            % Use sorted unique numbers
            cellIdx = unique(cellNumbers);
        end

        % cellIdx = [14];
        for i = 1:numel(cellIdx)
            %% Get Donor Excitation Donor Emission (DD) images
            DDIdx = find(contains({DDFiles.name}, [sampleName num2str(i) separator]));
            DDName = DDFiles(DDIdx).name;
            data.DD = imread([dataPath filesep DDName]);

            %% Get Donor Excitation Acceptor Emission (DA) images
            DAIdx = find(contains({DAFiles.name}, [sampleName num2str(i) separator]));
            DAName = DAFiles(DAIdx).name;
            data.DA = imread([dataPath filesep DAName]);

            %% Get Acceptor Excitation Donor Emission (AD) images
            ADIdx = find(contains({ADFiles.name}, [sampleName num2str(i) separator]));
            ADName = ADFiles(ADIdx).name;
            data.AD = imread([dataPath filesep ADName]);

            %% Get Acceptor Excitation Acceptor Emission (AA) images
            AAIdx = find(contains({AAFiles.name}, [sampleName num2str(i) separator]));
            AAName = AAFiles(AAIdx).name;
            data.AA = imread([dataPath filesep AAName]);

            %% Get Ilastik mask
            % Get list of possible files
            searchStr = [sampleName num2str(i) separator '*_FASegment.tiff'];
            disp('Constructed mask search string:')
            disp(searchStr)
            
            % Find matching mask file
            maskFiles = dir(fullfile(maskPath, segmentfolder, searchStr)); % Find matching mask file
            maskName = maskFiles.name;
            fprintf('Using mask: %s\n', maskName);
            data.mask = imread(fullfile(maskPath, segmentfolder, maskName));
                
            %% get name
            % str = [separator DDTag];
            % idx = strfind(DDName, str);
            % DDName(idx:idx+length(str)-1) = [];
            % 
            % str = dataType;
            % idx = strfind(DDName, str);
            % DDName(idx:idx+length(str)-1) = [];    
            data.name = DDName;

            % parse
            data.dataPath = dataPath;
            data.outputPath = outputPath;
            data.sampleIndex = [sampleName num2str(i)];
            data.index = i;
            %data.logicalSize = 700; 

            %% ROIs
            if contains(SampleFolder{E}, 'Vin')
                minSize = user.minSize_Vin;
            else
                minSize = user.minSize_VE;
            end
            data.mask = ROIandSeg(data, flags, minSize); % data.mask
            
            %% FRETint analysis
            close all;
            outputMetrics = IntFRET(data); % data.mask, data.DD etc.

            % Add identifiers
            outputMetrics.CellID     = i;
            % Append
            allOutputMetrics(i,:) = outputMetrics;
            save(fullfile(outputPath, 'allOutputMetrics.mat'), 'allOutputMetrics');
        end
        %close all;
    end
end   
            
disp('FRETint completed')
close all;





%% OLD VERSION
% %renameFoldersorFiles(rootDir, targetType, oldName, newName, lookName, filterName)
% %renameFoldersorFiles('D:\3_FRETTS_HBMVEC_analysis\2026_03_2222_FRETint', 'file', 'RAW_ch00', 'DD_ch00', '', 'RAW');
% 
% % -------------------------------------------------------------
% % Bulk FRET processing: masks in SegmentD_4, channels in parent
% % Cell discovery handles names with "CellX" or bare "X" (numbers).
% % -------------------------------------------------------------
% clear Tall;
% FolderInput   = 'D:\3_FRETTS_HBMVEC_analysis\2026_03_2222_FRETint\20260322_VinTS_FL_CT_si72';
% MaskSubfolder = fullfile(FolderInput, 'SegmentD_4');
% 
% if ~isfolder(FolderInput)
%     error('FolderInput does not exist: %s', FolderInput);
% end
% if ~isfolder(MaskSubfolder)
%     error('Mask subfolder not found: %s', MaskSubfolder);
% end
% 
% % -------------------------------------------------------------
% % Discover cells by scanning ALL files in SegmentD_4
% % We try to extract a cell number via either:
% %   1) "Cell123"    --> 123
% %   2) a standalone number "123" (not part of a longer digit sequence)
% % -------------------------------------------------------------
% maskFiles = dir(fullfile(MaskSubfolder, '*'));
% maskFiles = maskFiles(~[maskFiles.isdir]);  % files only
% 
% extractCellNum = @(name) local_extract_cell_number(name);
% 
% cellNums = [];
% for k = 1:numel(maskFiles)
%     n = extractCellNum(maskFiles(k).name);
%     if ~isnan(n)
%         cellNums(end+1) = n;
%     end
% end
% 
% cellNums = unique(cellNums);
% 
% if isempty(cellNums)
%     error('No cell numbers could be extracted from files in %s.', MaskSubfolder);
% end
% 
% fprintf('Found %d unique cells (from masks in SegmentD_4).\n', numel(cellNums));
% 
% % -------------------------------------------------------------
% % Prepare master table
% % -------------------------------------------------------------
% Tall = table();
% 
% % -------------------------------------------------------------
% % Process each cell
% % -------------------------------------------------------------
% for c = 1:numel(cellNums)
%     cellNum = cellNums(c);
%     cellID  = sprintf('Cell%d', cellNum);
%     fprintf('Processing %s ...\n', cellID);
% 
%     % -------------------------------
%     % Find the mask file for this cell number (first match)
%     % Strategy: match if filename contains "Cell<Num>" OR bare "<Num>"
%     % -------------------------------
%     maskPath = local_find_first_for_cell(MaskSubfolder, cellNum, {'FASegment'}, true);
%     if isempty(maskPath)
%         warning('No mask found for %s in SegmentD_4. Skipping.', cellID);
%         continue;
%     end
% 
%     % -------------------------------
%     % Find ch1 (ch00) and ch2 (ch01) in FolderInput
%     % Search allowing both "Cell<Num>" and bare "<Num>" patterns
%     % -------------------------------
%     ch1Path = local_find_first_for_cell(FolderInput, cellNum, {'ch00'}, false);
%     ch2Path = local_find_first_for_cell(FolderInput, cellNum, {'ch01'}, false);
% 
%     if isempty(ch1Path) || isempty(ch2Path)
%         warning('Missing channels for %s (ch00 or ch01). Skipping.', cellID);
%         continue;
%     end
% 
%     % -------------------------------
%     % Load images
%     % -------------------------------
%     try
%         mask = imread(maskPath) > 0;  % ensure binary
%     catch ME
%         warning('Failed to read mask for %s: %s. Skipping.', cellID, ME.message);
%         continue;
%     end
% 
%     try
%         ch1 = imread(ch1Path);
%         ch2 = imread(ch2Path);
%     catch ME
%         warning('Failed to read channels for %s: %s. Skipping.', cellID, ME.message);
%         continue;
%     end
% 
%     % Defensive: ensure sizes match
%     if ~isequal(size(mask), size(ch1), size(ch2))
%         warning('%s: Size mismatch (mask: %s, ch1: %s, ch2: %s). Skipping.', ...
%             cellID, mat2str(size(mask)), mat2str(size(ch1)), mat2str(size(ch2)));
%         continue;
%     end
% 
%     % -------------------------------
%     % Run your computeFRET (assumed on path)
%     % -------------------------------
%     T = IntFRET(mask, ch1, ch2); % MAIN FUNCTION
% 
%     % Add identifiers
%     T.CellID     = {cellID};
%     T.CellNumber = cellNum;
% 
%     % Append
%     Tall = [Tall; T];
% 
%     % % Access the FRET image:
%     % F = T.FRETimage{c};
%     % %figure; imagesc(F); axis image off; colorbar; title('FRET ratio (ch2/ch1)'); xlim([0 user.FRET_ColorbarLim]);
%     % figure;
%     % imagesc(F, [0 1]);   % Force the displayed range to 0–1
%     % axis image off;
%     % colormap(jet);       % Or 'turbo', 'parula', etc.
%     % colorbar;
%     % caxis([0 1]);        % Ensures the colorbar also spans 0–1
%     % title('FRET ratio (0–1 scaled display)');
% end
% 
% % -------------------------------------------------------------
% % Sort by CellNumber (optional)
% % -------------------------------------------------------------
% if ~isempty(Tall)
%     [~, ord] = sort(Tall.CellNumber);
%     Tall = Tall(ord, :);
% end
% 
% fprintf('Processing complete. %d cells added to table.\n', height(Tall));
% 
% % Optional: save results
% % writetable(Tall, fullfile(FolderInput, 'FRET_results.xlsx'));  % tables only (no arrays)
% % save(fullfile(FolderInput, 'FRET_results.mat'), 'Tall');        % preserves arrays in table
% 
% % ------------------------ Local helpers -----------------------
% function n = local_extract_cell_number(fname)
% % Extract a cell number from filename using robust patterns.
% % Priority:
% %   1) "Cell<digits>"
% %   2) standalone "<digits>" not immediately adjacent to other digits
%     n = NaN;
% 
%     % Pattern 1: "Cell123"
%     tok = regexp(fname, 'Cell(\d+)', 'tokens', 'once');
%     if ~isempty(tok)
%         n = str2double(tok{1});
%         return;
%     end
% 
%     % Pattern 2: standalone number (not part of a longer digit run)
%     % We use lookarounds to ensure non-digit boundaries
%     toks = regexp(fname, '(?<!\d)(\d+)(?!\d)', 'tokens');
%     if ~isempty(toks)
%         % If multiple numbers, pick the first; adjust if you prefer last
%         n = str2double(toks{1}{1});
%     end
% end
% 
% function fpath = local_find_first_for_cell(folder, cellNum, mustContainList, restrictToMasks)
% % Find the first file path in 'folder' that appears to belong to cellNum.
% % Matches either "Cell<Num>" OR a standalone "<Num>" in the filename.
% % Additionally, filename must contain all substrings in mustContainList.
% % If restrictToMasks = true, we prioritize files that look like masks (contain 'Segment' or 'FASegment').
% 
%     fpath = '';
%     files = dir(folder);
%     files = files(~[files.isdir]);
% 
%     % Precompile patterns
%     patCell  = sprintf('Cell%d', cellNum);
%     patBare  = sprintf('(?<!\\d)%d(?!\\d)', cellNum);  % standalone number
% 
%     for pass = 1:2
%         % pass 1: require "Cell<Num>"
%         % pass 2: allow standalone "<Num>"
%         for k = 1:numel(files)
%             name = files(k).name;
% 
%             % Must contain all required substrings (like 'ch00' or 'FASegment')
%             if ~isempty(mustContainList)
%                 ok = all(cellfun(@(s) contains(name, s), mustContainList));
%                 if ~ok, continue; end
%             end
% 
%             % If restricting to mask-like names, enforce a hint
%             if restrictToMasks
%                 if ~(contains(name, 'Segment') || contains(name, 'FASegment'))
%                     continue;
%                 end
%             end
% 
%             switch pass
%                 case 1
%                     if ~contains(name, patCell), continue; end
%                 case 2
%                     if isempty(regexp(name, patBare, 'once')), continue; end
%             end
% 
%             fpath = fullfile(files(k).folder, name);
%             return;
%         end
%     end
% end
