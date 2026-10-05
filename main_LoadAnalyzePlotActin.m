%% Load Actin staining .tif files (ch01) into table T_Act
clear; close all; clc;

% ====================================================================
%% ================= 0_User Configuration =================
% ====================================================================
Input.mainDir = 'D:\3_FRETTS_HBMVEC_analysis\2026_07_23_ActinStaining';

% -- Which folders to read? ({} = read everything found automatically) --
Input.ConditionFolders = {};   % e.g. {'CT_1kPa', 'CCM_10kPa'}
Input.SubFolders       = {};   % e.g. {'2026_03_303101_VE', '2026_04_11_Vin'} (applies to all condition folders)

% Optional: different subfolder selection for specific condition folders.
% Overrides Input.SubFolders for that condition folder only.
Input.SubFoldersPerCondition = containers.Map('KeyType', 'char', 'ValueType', 'any');
% Input.SubFoldersPerCondition('CT_1kPa')  = {'2026_03_303101_VE'};
% Input.SubFoldersPerCondition('CCM_10kPa') = {'2026_04_11_Vin', '2026_04_11_VE'};

% -- File selection --
Input.dataType   = '.tif';
Input.channelTag = 'ch01';     % only files containing this tag are read (actin)
Input.ImageSize  = 512;        % expected image size (only used for a warning)

% -- Category orders (ordinal categoricals, same idea as Order_TableByRules) --
Input.proteinOrder   = {'Vin', 'VECad'};
% Folder-name suffix -> protein name used in the table (folders end in _VE, table says VECad)
Input.proteinMap = containers.Map({'vin', 've', 'vecad'}, {'Vin', 'VECad', 'VECad'});
Input.diseaseOrder   = {'CT', 'CCM'};
Input.extrainfoOrder = {'1kPa', '10kPa'};
Input.sampletypeOrder = {'TS', 'TL'};   % TS = tension sensor, TL = tension-less control; defaults to TS if not in folder name

% ====================================================================
%% ================= 1_LOAD DATA =================
% ====================================================================
% --- Condition folders (e.g. CT_1kPa) ---
if isempty(Input.ConditionFolders)
    ConditionFolders = listSubfolders(Input.mainDir);
else
    ConditionFolders = Input.ConditionFolders;
end

% Pre-allocate containers (one entry per image, grown as we go)
DataFolder = strings(0,1);
SampleFolder = strings(0,1);
sampletype = strings(0,1);
disease = strings(0,1);
extrainfo = strings(0,1);
dateStr = strings(0,1);
sampleTypeStr = strings(0,1);
FileName = strings(0,1);
ActinImage = cell(0,1);
nRGB = 0; nRGBdiff = 0;   % counters for RGB-stored tifs
ImgMin = zeros(0,1); ImgMax = zeros(0,1); SatFrac = zeros(0,1);   % intensity QC per image
for c = 1:numel(ConditionFolders)
    condName = ConditionFolders{c};
    condPath = fullfile(Input.mainDir, condName);
    if ~isfolder(condPath)
        warning('Condition folder not found, skipping: %s', condPath);
        continue
    end

    % Parse e.g. 'CT_1kPa' -> disease 'CT', extrainfo '1kPa'
    tokCond = regexp(condName, '^([^_]+)_(.+)$', 'tokens', 'once');
    if isempty(tokCond)
        warning('Cannot parse disease/extrainfo from "%s", skipping.', condName);
        continue
    end
    thisDisease   = string(tokCond{1});
    thisExtrainfo = string(tokCond{2});

    % --- Subfolders (e.g. 2026_03_303101_VE) ---
    if isKey(Input.SubFoldersPerCondition, condName)
        SubFolders = Input.SubFoldersPerCondition(condName);
    elseif ~isempty(Input.SubFolders)
        SubFolders = Input.SubFolders;
    else
        SubFolders = listSubfolders(condPath);
    end

    for s = 1:numel(SubFolders)
        subName = SubFolders{s};
        subPath = fullfile(condPath, subName);
        if ~isfolder(subPath)
            warning('Subfolder not found, skipping: %s', subPath);
            continue
        end

        % Parse e.g. '2026_03_303101_VE' -> date '2026_03_303101', protein 'VE'
        tokSub = regexp(subName, '^(.+)_(VECad|VE|Vin)$', 'tokens', 'once', 'ignorecase');
        if isempty(tokSub)
            warning('Cannot parse date/protein from "%s", skipping.', subName);
            continue
        end
        thisDate = string(tokSub{1});
        % map folder suffix (VE / VECad / Vin) to the canonical name (VECad / Vin)
        thisProtein = string(Input.proteinMap(lower(tokSub{2})));

        % TS / TL from the folder names (e.g. '..._VinTS', '..._TL_...'); default TS if absent
        tsTl = regexp([condName '_' subName], '(?<![A-Za-z])(?:Vin|VECad|VE)?(TS|TL)(?![A-Za-z])', 'tokens');
        if isempty(tsTl)
            thisSampleType = "TS";
        else
            found = unique(cellfun(@(c) c{1}, tsTl, 'UniformOutput', false), 'stable');
            if numel(found) > 1
                warning('Both TS and TL found in "%s/%s"; using %s.', condName, subName, found{1});
            end
            thisSampleType = string(found{1});
        end

        % --- Files: only actin (ch01) .tif ---
        files = dir(fullfile(subPath, ['*' Input.channelTag '*' Input.dataType]));
        if isempty(files)
            warning('No "%s" files in %s', Input.channelTag, subPath);
            continue
        end

        for f = 1:numel(files)
            filePath = fullfile(files(f).folder, files(f).name);
            img = imread(filePath);
            if ndims(img) > 2
                % RGB (or RGBA) tif: keep only the plane(s) that contain signal
                img = img(:,:,1:min(3,size(img,3)));
                nRGB = nRGB + 1;
                activePlanes = find(squeeze(max(max(img,[],1),[],2)) > 0);
                if isempty(activePlanes)
                    img = img(:,:,1);                 % empty image, keep first plane
                elseif isscalar(activePlanes)
                    img = img(:,:,activePlanes);      % pseudo-coloured single channel (e.g. only red)
                elseif isequal(img(:,:,1), img(:,:,2)) && isequal(img(:,:,1), img(:,:,3))
                    img = img(:,:,1);                 % grayscale stored as RGB
                else
                    img = max(img, [], 3);            % signal in several planes
                    nRGBdiff = nRGBdiff + 1;
                    warning('%s: signal in several RGB planes; using max over planes.', files(f).name);
                end
            end
            % --- intensity QC on the native-class image (before double conversion) ---
            satValue = double(intmax(class(img)));          % 255 for uint8, 65535 for uint16
            thisSat  = mean(img(:) == satValue);            % fraction of saturated pixels
            ImgMin(end+1,1)  = double(min(img(:)));
            ImgMax(end+1,1)  = double(max(img(:)));
            SatFrac(end+1,1) = thisSat;
            if thisSat > 0.001
                warning('%s: %.2f%% saturated pixels.', files(f).name, 100*thisSat);
            end

            img = double(img);
            if ~isequal(size(img), [Input.ImageSize Input.ImageSize])
                warning('%s is %dx%d, expected %dx%d.', files(f).name, ...
                    size(img,1), size(img,2), Input.ImageSize, Input.ImageSize);
            end

            DataFolder(end+1,1)   = string(subName);
            SampleFolder(end+1,1) = string(condName);
            sampletype(end+1,1)   = thisProtein;
            disease(end+1,1)      = thisDisease;
            extrainfo(end+1,1)    = thisExtrainfo;
            dateStr(end+1,1)      = thisDate;
            sampleTypeStr(end+1,1) = thisSampleType;
            FileName(end+1,1)     = string(files(f).name);
            ActinImage{end+1,1}   = img;
        end
        fprintf('Loaded %2d images: %s / %s\n', numel(files), condName, subName);
    end
end

% ====================================================================
%% ================= 2_BUILD TABLE =================
% ====================================================================
if nRGB > 0
    fprintf('%d of %d images were stored as RGB; %d had signal in several planes (max over planes used).\n', ...
        nRGB, numel(ActinImage), nRGBdiff);
end
dateLevels = unique(dateStr);   % sorted -> chronological for YYYY_MM_... names

T_Act = table(DataFolder, SampleFolder, sampletype, disease, extrainfo, ActinImage, ...
    'VariableNames', {'DataFolder','SampleFolder','sampletype','disease','extrainfo','ActinImage'});

T_Act.dateCat      = categorical(dateStr,   cellstr(dateLevels),         'Ordinal', true);
T_Act.proteinCat   = categorical(sampletype, Input.proteinOrder,   'Ordinal', true);
T_Act.diseaseCat   = categorical(disease,    Input.diseaseOrder,   'Ordinal', true);
T_Act.extrainfoCat = categorical(extrainfo,  Input.extrainfoOrder, 'Ordinal', true);
T_Act.sampletypeCat = categorical(sampleTypeStr, Input.sampletypeOrder, 'Ordinal', true);   % TS / TL (needed by plotMetricScatter)
T_Act.expCat       = T_Act.dateCat;   % biological replicate = experiment date (used by plotMetricScatter / stats)
T_Act.FileName     = FileName;   % extra: handy for tracing back to the raw file
T_Act.ImgMin       = ImgMin;     % QC: min / max raw value and fraction of saturated pixels
T_Act.ImgMax       = ImgMax;
T_Act.SatFrac      = SatFrac;
fprintf('Saturation: %d of %d images have >0.1%% saturated pixels (max SatFrac = %.2f%%).\n', ...
    sum(SatFrac > 0.001), numel(SatFrac), 100*max(SatFrac));

% Sanity check: any value that did not match the expected orders becomes <undefined>
if any(ismissing(T_Act(:, {'proteinCat','diseaseCat','extrainfoCat','sampletypeCat'})), 'all')
    warning('Some rows have <undefined> categories; check the Input.*Order lists.');
end

% Overview: number of images per combination
disp(groupsummary(T_Act, {'dateCat','proteinCat','diseaseCat','extrainfoCat'}));

%% Save (optional)
save('C:\Users\sarah\OneDrive - KU Leuven\Analysis_Codes_Matlab\Code saveguarding Jan26\20260114_FDA_V3_SarahBulk\6_SavedData\TAct.mat', 'T_Act', '-v7.3')

% ====================================================================
%% ================= 2_RUN ANALYSIS =================
% ====================================================================
T_Act = Calc_ActinIntensity(T_Act, Threshold=20);       % fixed threshold
%T_Act = Calc_ActinIntensity(T_Act, Method="isodata");  % automatic, threshold per image

% ====================================================================
%% ================= 3_PLOT =================
% ====================================================================
% Input
Input.ContrastExcludePatterns = {{'kPa', 'R1'},{'kPa', 'R2'},{'TL', 'TS'}}; % comparisons excluded from stat because they make no sense
Input.ContrastRemove = {{'CT_1kPa_TS','CCM_10kPa_TS'}, {'CT_10kPa_TS','CCM_1kPa_TS'}};
Input.condColors = containers.Map(...
    {'CT_rest',  'CT_R1',    'CT_R2',    'CT_10kPa',    'CT_1kPa', ...
     'CCM_rest', 'CCM_R1',   'CCM_R2',   'CCM_10kPa',   'CCM_1kPa'}, ...
    {'blues',     'indigos',   'indigos',   'bluegreens',   'bluegreens', ...
     'redes',     'magentas',  'magentas',  'redoranges',   'redoranges'});
Input.CondOrderGel = {'CT_1kPa_TS', 'CT_10kPa_TS', 'CCM_1kPa_TS', 'CCM_10kPa_TS'};
Input.OutPath = 'C:\Users\sarah\OneDrive - KU Leuven\WRITING_Articles\PaperFRETTS-CCM\1_FigureFiles';

%% Plot
plotMetricScatter(T_Act, 'AverageIntensity', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors,...
    ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_Actin_AvInt'), 'svg', Width=16, Height=12);

%% ====================================================================
% Local functions
% ====================================================================
function names = listSubfolders(parentPath)
% Names of all subfolders in parentPath (without '.' and '..')
    d = dir(parentPath);
    d = d([d.isdir]);
    names = {d.name};
    names = names(~ismember(names, {'.', '..'}));
end