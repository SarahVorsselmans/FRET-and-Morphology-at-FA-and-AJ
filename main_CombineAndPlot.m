%% Combine all allOutputMetrics.mat files, analyze & make figures
clear; close all; clc;

%% When DIP Image is really messing everything up:
% restoredefaultpath
% rehash toolboxcache

% ====================================================================
%% ================= 0_User Configuration =================
% ====================================================================
% -- Image inputs --
Input.Image_Resolution = 512; % 512x512 image size
Input.Image_ImSize = 123.2577; % µm
Input.PixSize = 0.24073767; % µm
Input.MinPixelFA = 5; % µm

% -- Categories & metrics --
Input.FACatCategories = {
    {'MinMajRatio', '<',  0.3},     ... % Cat1
    {'MinMajRatio', '<',  0.5},     ... % Cat2
    {'MinMajRatio', '>=', 0.5, 'Area', '<',  25},   ... % Cat3
    {'MinMajRatio', '>=', 0.5, 'Area', '>=', 25}};      % Cat4
% Input.FACatCategories2 = {
%     {'MinMajRatio', '<',  0.25, 'MajAx', '>=',  8},    ... % linear - 1
%     {'MinMajRatio', '<',  0.5,  'MajAx', '>=',  8},     ... % ellipsoid - 2
%     {'MajAx', '<',  8},     ... % small - 3
%     {'MinMajRatio', '>=', 0.5}};        % big round - 4
    
Input.AJCatCategories = {
    {'Thickness','>=',4, 'EulerNum','<',0}  ... % Cat1 – thick & holey
    {}};     % Cat2 – thin & uniform % 'Thickness','<',5,'EulerNum','==',1

Input.CatNames = {'Linear', 'Elliptical', 'Small Circular', 'Large Circular'};
Input.CatNamesAJ = {'Reticular', 'Linear'};
% Input.CondOrder = {'CT_rest', 'CT_R1', 'CT_R2','CT_10kPa', 'CT_1kPa','CCM_rest', 'CCM_R1', 'CCM_R2','CCM_10kPa', 'CCM_1kPa'};
% Input.CondOrder2 = {'CT_rest', 'CCM_rest', 'CT_R1', 'CCM_R1', 'CT_R2', 'CCM_R2', 'CT_10kPa', 'CCM_10kPa', 'CT_1kPa', 'CCM_1kPa'};
Input.CondOrder = {'CT_rest_TS', 'CT_R1_TS', 'CT_R2_TS','CT_10kPa_TS', 'CT_1kPa_TS','CCM_rest_TS', 'CCM_R1_TS', 'CCM_R2_TS','CCM_10kPa_TS', 'CCM_1kPa_TS', 'CT_rest_TL', 'CCM_rest_TL'};
Input.CondOrder2 = {'CT_rest_TS', 'CCM_rest_TS', 'CT_R1_TS', 'CCM_R1_TS', 'CT_R2_TS', 'CCM_R2_TS', 'CT_10kPa_TS', 'CCM_10kPa_TS', 'CT_1kPa_TS', 'CCM_1kPa_TS', 'CT_rest_TL', 'CCM_rest_TL'};
Input.SpringOrder = {'TS_F40', 'TL_F40', 'TS_FL', 'TL_FL','TS_stHP35', 'TL_stHP35'};
Input.CondOrderGel = {'CT_1kPa_TS', 'CT_10kPa_TS', 'CCM_1kPa_TS', 'CCM_10kPa_TS'};

Input.condColors = containers.Map(...
    {'CT_rest',  'CT_R1',    'CT_R2',    'CT_10kPa',    'CT_1kPa', ...
     'CCM_rest', 'CCM_R1',   'CCM_R2',   'CCM_10kPa',   'CCM_1kPa'}, ...
    {'blues',     'indigos',   'indigos',   'bluegreens',   'bluegreens', ...
     'redes',     'magentas',  'magentas',  'redoranges',   'redoranges'});
Input.condColorsClass = containers.Map(...
    Input.CatNames, ...
    {'bluepurple_custom', 'green_custom', 'yellow_custom', 'purple_custom'});
Input.condColorsClassAJ = containers.Map(...
    Input.CatNamesAJ, ...
    {'bluepurple_custom', 'green_custom'});

Input.CompareVariablesFA = {'FA_Count', 'Av_AlignIndex', 'FracCat1', 'FracCat4','Av_AverageIntensity', 'FRETav'};
Input.CompareVariablesAJ = {'CellAv_MinMajRatio', 'Thickness_Skeleton_Av', 'FracCat1', 'FracCat2', 'Av_AverageIntensity', 'FRETav'};

Input.CompareVariables = {'FA_FRETav', 'FA_FA_Count', 'FA_Av_PixSize','FA_Av_MinMajRatio',... 
    'FA_Av_MajAx', 'FA_Av_FADens_r5', 'FA_Av_AverageIntensity', 'FA_Av_AlignIndex', 'FA_FracCat1', 'FA_FracCat4',...
    'FRETav', 'CellAvPixSize', 'CellAv_MajAx', 'CellAv_MinMajRatio', 'CellAlignIndex', ...
    'Thickness_Skeleton_Av', 'Thickness_Area_Av', 'Av_AverageIntensity', 'FracCat1', 'FracCat2'};

% -- Stat --
Input.ContrastExcludePatterns = {{'kPa', 'R1'},{'kPa', 'R2'},{'TL', 'TS'}, {'rest', 'kPa'}}; % comparisons excluded from stat because they make no sense
% {'TL', 'kPa'},{'TL', 'R1'},{'TL', 'R2'}
Input.ContrastRemove = {{'CT_1kPa_TS','CCM_10kPa_TS'}, {'CT_10kPa_TS','CCM_1kPa_TS'}};

% -- Path inputs --
Input.mainDir = 'D:\3_FRETTS_HBMVEC_analysis';
Input.OutPath = 'C:\Users\sarah\OneDrive - KU Leuven\WRITING_Articles\PaperFRETTS-CCM\1_FigureFiles';
%Input.mainDir = 'C:\Users\sarah\OneDrive - KU Leuven\Analysis_Codes_Matlab\TEMPORARY_VECadtests\';
Input.DataFolders = {
    % '2025_08_11_StableVinCCM_fix72h','2025_08_18_StableVECadVin_fix72h',...
    % '2025_09_01_StableVinVE_4h24hgel','2025_09_22_TLsAndVEFL',...
    % '2025_10_0506_VinVETSM_GelR12', '2025_10_1314_VinVE_TSTL_R1R2_48h7d', ...
    % '2025_10_2627_VinVE_GelROCKLive', '2025_11_0206_VinVEGel_siROCK', ...
    % '2025_11_1011_VinVETSTL', '2025_11_1616_VinVETS_FL_48h_live', ...
    % '2025_11_1718_VinVETS_FL_72h_fix', '2025_11_2424_VinVEGel', 
    % '2026_01_252627_WellRockGel','2026_02_310102_WellRockGel', ...
    % '2026_03_0809_WellRockGel', '2026_03_1616_VinVE24welltest',...
    % '2026_03_303101_WellRockGel', '2026_04_1115_OldTests',...
    % '2026_04_111415_WellRockGel', '2026_04_181920_WellRockGel',...
    % '2026_04_272727_VinVEGel_20260124_ok','2026_04_292929_VinVEGel_20260130_ok'
    % APPROVED:
    '2025_08_11_StableVinCCM_fix72h_ok',...
    '2025_08_18_StableVECadVin_fix72h_ok',...
    '2025_09_22_TLsAndVEFL_ok',...
    '2025_11_1718_VinVETS_FL_72h_fix_ok',...
    '2025_11_1616_VinVETS_FL_48h_live_ok',...
    '2026_03_303101_WellRockGel_ok',...
    '2026_04_111415_WellRockGel_ok',...
    '2026_04_181920_WellRockGel_ok',...
    '2026_04_272727_VinVEGel_20260124_ok',...
    '2026_04_292929_VinVEGel_20260130_ok'};

% Format SampleFolders: YYYYMMDD_<Protein><TS/TL>_<spring>_<disease>_<extrainfo>
% Metrics to extract (extensible)
Input.numFields   = {}; % fields to convert to numeric % should fix this in FRETcode at some point 
Input.cellFields  = {'FRETimageFA', 'IntA', 'IntD', 'FRETav'}; % fields to keep as-is (cell column)

% !! Additional
% Input.includePatterns_VinFL = {'VinTS_FL', 'VinTL_FL'}; % CAREFUL! : include with ... ! or ! with ...
% Input.excludePatterns_Vin = {'ish', 'si48', '20250922','20251116'};
% T_Vin = Filter_IncludeExclude(T_Vin, Input.includePatterns_VinFL, Input.excludePatterns_Vin);
% %
Input.includePatterns_VECadFL = {'VECadTS_FL', 'VECadTL_FL'};
Input.excludePatterns_VECad = {'ish', 'si48', '20250922', '20251116'}; % STILL NECESSARY!
T_VECad = Filter_IncludeExclude(T_VECad, Input.includePatterns_VECadFL, Input.excludePatterns_VECad);

% ====================================================================
%% ================= 1_LOAD DATA =================
% ====================================================================
allData = Load_DataFolders(Input, '');
%%
% save('D:\SarahV\20260216_AllData.mat', 'allData', '-v7.3')
save('C:\Users\sarah\OneDrive - KU Leuven\Analysis_Codes_Matlab\Code saveguarding Jan26\20260114_FDA_V3_SarahBulk\6_SavedData\20260810_TVECad.mat', 'T_VECad', '-v7.3')
% save('C:\Users\sarah\OneDrive - KU Leuven\Analysis_Codes_Matlab\Code saveguarding Jan26\20260114_FDA_V3_SarahBulk\6_SavedData\20260604_TVinSpring.mat', 'T_VinSpring', '-v7.3')
% save('C:\Users\sarah\OneDrive - KU Leuven\Analysis_Codes_Matlab\Code saveguarding Jan26\20260114_FDA_V3_SarahBulk\6_SavedData\20260604_TVinLiveFixed.mat', 'T_VinLiveFixed', '-v7.3')
% save('C:\Users\sarah\OneDrive - KU Leuven\Analysis_Codes_Matlab\Code saveguarding Jan26\20260114_FDA_V3_SarahBulk\6_SavedData\20260604_TVECad.mat', 'T_VECad', '-v7.3')
% save('C:\Users\sarah\OneDrive - KU Leuven\Analysis_Codes_Matlab\Code saveguarding Jan26\20260114_FDA_V3_SarahBulk\6_SavedData\20260604_TVECadSpring.mat', 'T_VECadSpring', '-v7.3')
% save('C:\Users\sarah\OneDrive - KU Leuven\Analysis_Codes_Matlab\Code saveguarding Jan26\20260114_FDA_V3_SarahBulk\6_SavedData\20260604_TVECadLiveFixed.mat', 'T_VECadLiveFixed', '-v7.3')

% ====================================================================
%% ================= 2_PRE-PROCESS and ORDER DATA =================
% ====================================================================
%% -- Remove additional regions --
% See Main_DataForElodie.m

%% -- Pre-process data --
% 1) Make Mask from FRETimageFA
for r = 1:height(allData)
    allData.Mask{r} = ~isnan(allData.FRETimageFA{r});
end

% 2) Mask Int image from intensity A 
% (+ make Unit16 into double and bg becomes NaN)
allData.Int = cellfun(@(m, a) ...
    setMaskedNaN(double(a) .* double(m)), ...
    allData.Mask, allData.IntA, ...
    'UniformOutput', false);
% just make into double: T_VECad.Int = cellfun(@double, T_VECad.Int, 'UniformOutput', false);
% Helper function
function y = setMaskedNaN(y)
    y(y == 0) = NaN;
end

% 3) Remove small regions from Mask for vinculin data
for i = 1:height(T_Vin)
    % Extra small-object removal, only for samples whose SampleFolder name contains 'Vin'
    if contains(T_Vin.SampleFolder{i}, 'Vin', 'IgnoreCase', true)
        T_Vin.Mask{i} = bwareaopen(T_Vin.Mask{i}, Input.MinPixelFA, 8); % removes areas smaller than X pixels
    end
end

%% -- Order & Categorize -- (creates ORDINAL categorical columns)
Input.springOrder     = {'F40','FL','stHP35'};
Input.proteinOrder    = {'Vin','VECad'};
Input.sampletypeOrder = {'TS','TL'};   % TS before TL
Input.diseaseOrder    = {'CT','CCM'};
Input.extrainfoOrder  = {'rest','R1','R2','10kPa','1kPa'};
Input.LiveFixedOrder  = {'Fixed','Live'};

allData_sorted = Order_TableByRules(allData, Input);
%T_Vin = Order_TableByRules(T_Vin, Input);
%T_VECad = Order_TableByRules(T_VECad, Input);
%T_VECadFL = Order_TableByRules(T_VECadFL, Input);
% T_VinLiveFixed   = Order_TableByRules(T_VinLiveFixed, Input);
% T_VECadLiveFixed = Order_TableByRules(T_VECadLiveFixed, Input);

%% -- Build ExportT from sorted T for GraphPad use --
% [ExportT] = ExportTGraphPad(allData, 'FRETav'); % adjust input

% ====================================================================
%% ================= 3_FILTER DATA =================
% ====================================================================
%% Vin
Input.includePatterns_VinFL = {'VinTS_FL', 'VinTL_FL'}; % CAREFUL! : include with ... ! or ! with ...
Input.excludePatterns_Vin = {'ish', 'si48', '20250922','20251116'};
T_Vin = Filter_IncludeExclude(allData_sorted, Input.includePatterns_VinFL, Input.excludePatterns_Vin);
%%
Input.includePatterns_VinSpring = {'20250811_VinT', '20250922_VinT'};
T_VinSpring = Filter_IncludeExclude(allData_sorted, Input.includePatterns_VinSpring, Input.excludePatterns_Vin); % ! should not be gel or ROCK!
Input.includePatterns_VinLiveFixed = {'20251116_VinT', '20251117_VinT', '20251118_VinT'};
T_VinLiveFixed = Filter_IncludeExclude(allData_sorted, Input.includePatterns_VinLiveFixed, {});

%% VECad
Input.includePatterns_VECadFL = {'VECadTS_FL', 'VECadTL_FL'};
Input.excludePatterns_VECad = {'ish', 'si48', '20250922', '20251116'};
T_VECad = Filter_IncludeExclude(allData_sorted, Input.includePatterns_VECadFL, Input.excludePatterns_VECad);
%%
Input.includePatterns_VECadSpring = {'20250818_VECadT', '20250922_VECadT'};
T_VECadSpring = Filter_IncludeExclude(allData_sorted, Input.includePatterns_VECadSpring, Input.excludePatterns_VECad); % ! should not be gel or ROCK!
Input.includePatterns_VECadLiveFixed = {'20251116_VECadT', '20251118_VECadT'};
T_VECadLiveFixed = Filter_IncludeExclude(allData_sorted, Input.includePatterns_VECadLiveFixed, {});

% ====================================================================
%% ================= 4_CALC VIN DATA =================
% ====================================================================
%% Shape, Density/NND, Intensity
T_Vin = Calc_FA_AllMetrics(T_Vin, Input, FRET=true); 
% Output: all Nx1 doubles per row
% Output: CC_Obj, FRET, Ecc, MajAx, MinAx, MinMajRatio, Area, Orientation
% Output: centroid, NND, NND(3), FA density (radius 5µm) (in µm)
% Output: IntensityPerPixel, TotalInt, AverageIntensity,
% AverageIntensityNorm, MaxInt, MinInt 

%% Averaging metrics (only if min. 5 FAs per FOV)
T_Vin = Calc_FA_AveragingMetrics(T_Vin, Input, FRET=true); 
% Output: FRETav, FA_Count, Av_PixSize, Av_UmSize, Median_FA_Size_pix,
% Total_FA_Area_pix, Total_FA_Area_um, Av_MinMajRatio, Av_Ecc, Av_NND3, 
% Av_FADens_r5, Av_MajAx, Av_MajAx_um, Av_AverageIntensity, 
% Av_MeanAngle, Av_CircVar, Av_AlignIndex

%% Ratio FRET & intensity
T_Vin = Calc_FA_TensionPerDensity(T_Vin);

%% FA categories
T_Vin = Calc_SplitFAsByShapeCategory(T_Vin, Input.FACatCategories, FRET=true); % FRET=true
% Output: MaskClass (image), Int_ & FRET_ for 'Cat1, Cat2, Cat3...'
%T_Vin2 = Calc_SplitFAsByShapeCategory(T_Vin, Input.FACatCategories2, FRET=true);

%% Divide table by categories and calculate metrics & averaging metrics
%CatTables = CalcSplit_FullTableByCategory(T_Vin, Input, FRET=true); % FRET=true !!!!
% Output: Structure with tables per category

%% Calculate Cohen's d between conditions for heatmap
%allcomparisonResultsFA = Calc_compareConditionsToReference_ALL(CatTables, Input.CompareVariablesFA);

% Extract some categories first as individual columns:
T_Vin.FracCat1 = T_Vin.FA_summaryFrac(:, 1);
T_Vin.FracCat4 = T_Vin.FA_summaryFrac(:, 4);
%%
comparisonResultsFA = Calc_compareConditionsToReference(T_Vin, Input.CompareVariablesFA);
%%%
% T_diffVinkPa = Calc_DifferenceTable(T_Vin, {'FA_Count','FRETav', 'Av_AverageIntensity'}, '1kPa', '10kPa', ...
%     compareCat='diseaseCat', XCat='extrainfoCat', BioRepCol='expCat', excludePatterns={'TL'});
% T_diffVinkPa_Rev = Calc_DifferenceTable(T_Vin, {'FA_Count','FRETav', 'Av_AverageIntensity'}, '10kPa', '1kPa', ...
%     compareCat='diseaseCat', XCat='extrainfoCat', BioRepCol='expCat', excludePatterns={'TL'});
% T_diffVinkPa2 = Calc_DifferenceTable(T_Vin, {'FA_Count','FRETav', 'Av_AverageIntensity'}, 'CT', 'CCM', ...
%     compareCat='extrainfoCat', XCat='diseaseCat', BioRepCol='expCat', includePatterns={'kPa'}, excludePatterns={'TL'});

%%
keyboard;
%% Same for springs
% T_VinSpring = Calc_FA_AllMetrics(T_VinSpring, Input, FRET=true); 
% T_VinSpring = Calc_FA_AveragingMetrics(T_VinSpring, Input, FRET=true); 
% T_VinSpring = Calc_SplitFAsByShapeCategory(T_VinSpring, Input.FACatCategories, FRET=true);

%% Other analysis possible:
% - A. Size–intensity scaling: Test whether larger FAs have proportionally more intensity (constant average), superlinear, or sublinear scaling.
% - B. Rim vs. core enrichment: Compute rim/core mean intensity ratio per FA (already included as EdgeCoreRatio in the function). Compare distributions across conditions.
% - C. Shape–intensity coupling: Correlate AverageIntensity or IntensityPerPixel with MajorAxisLength, MinorAxisLength, and aspect ratio (Major/Minor). This tells whether elongated adhesions tend to be brighter.


% ====================================================================
%% ================= 5_CALC VECad DATA =================
% ====================================================================
%% CellMask
T_VECad = LOOPcomputeCellSize(T_VECad); % not the best function yet, lot of redundancy
%% Cell metrics
T_VECad = Calc_CellShapeMetrics(T_VECad);
T_VECad = Calc_Cell_AveragingMetrics(T_VECad, Input);

%% Thickness
[T_VECad] = Calc_Thickness(T_VECad, Input);

%% Split TS from TL
% Input.includePatterns_VECadTS = {'VECadTS_FL'};
% T_VECad = Filter_IncludeExclude(T_VECad, Input.includePatterns_VECadTS, Input.excludePatterns_VECad);

%% Classify VECad - paused for now
T_VECad = Calc_classifyAdhesions(T_VECad, AJCatCategories=Input.AJCatCategories);
T_VECad.FracCat1 = T_VECad.FA_summaryFrac(:, 1);
T_VECad.FracCat2 = T_VECad.FA_summaryFrac(:, 2);
%% Calculate FRET & Int per classification
T_VECad = Calc_AdhesionSkelCategoryStats(T_VECad, 'FRET', true);

%% Divide table by categories and calculate metrics & averaging metrics
% CatTablesAJ = CalcSplit_FullTableByCategory_AJ(T_VECad, Input); % FRET=true !!!!
% Output: Structure with tables per category

%% Calculate Cohen's d between conditions for heatmap
% allcomparisonResultsAJ = Calc_compareConditionsToReference_ALL(CatTablesAJ, Input.CompareVariablesAJ);
comparisonResultsAJ = Calc_compareConditionsToReference(T_VECad, Input.CompareVariablesAJ);
%%
T_diffVEkPa = Calc_DifferenceTable(T_VECad, {'Thickness_Skeleton_Av','FRETav', 'Av_AverageIntensity'}, '1kPa', '10kPa', ...
    compareCat='diseaseCat', XCat='extrainfoCat', BioRepCol='expCat', excludePatterns={'TL'});
T_diffVEkPa_Rev = Calc_DifferenceTable(T_VECad, {'Thickness_Skeleton_Av','FRETav', 'Av_AverageIntensity'}, '10kPa', '1kPa', ...
    compareCat='diseaseCat', XCat='extrainfoCat', BioRepCol='expCat', excludePatterns={'TL'});
T_diffVEkPa2 = Calc_DifferenceTable(T_VECad, {'Thickness_Skeleton_Av','FRETav', 'Av_AverageIntensity'}, 'CT', 'CCM', ...
    compareCat='extrainfoCat', XCat='diseaseCat', BioRepCol='expCat', includePatterns={'kPa'}, excludePatterns={'TL'});

%% Boris' PIC pore code
%[T_VECad] = computeThickness_PoreCode(T_VECad);
% if ismember('FRETindex', T.Properties.VariableNames)
%     % Preallocate: one cell per row for the diameter vectors
%     T.ThicknessBub   = cell(height(T), 1);
%     % Preallocate: one scalar per row for the averages
%     T.ThicknessBubAv = nan(height(T), 1);
% 
%     for i = 1:height(T)
%         % Guard against missing/empty masks
%         Mi = T.FRETindex{i};
%         % Replace invalid values with background (100)
%         Mi(~isfinite(Mi)) = 100;
%         % Build binary mask: 100 -> 0 (background), anything else -> 1 (foreground)
%         Mi = Mi ~= 100;
%         if isempty(Mi)
%             T.ThicknessBub{i}   = [];      % keep empty
%             T.ThicknessBubAv(i) = NaN;     % average is NaN
%             continue
%         end
% 
%         % Ensure logical for bwmorph (it expects a binary/logical image)
%         Mi = logical(Mi);
% 
%         % Run your function on the closed mask
%         pores2D = getPoreProps2D(bwmorph(Mi, 'close'));
% 
%         % Get the diameter vector (Nx1 double). Make sure it's a column.
%         diam = pores2D.diameter(:);
% 
%         % Store the vector into a cell and the average into a numeric column
%         T.ThicknessBub{i}   = diam;
%         T.ThicknessBubAv(i) = mean(diam, 'omitnan');
%     end
% end
keyboard;

%% Extra stats
% combo = string(T_Vin.diseaseCat) + "_" + string(T_Vin.extrainfoCat) + "_" + string(T_Vin.sampletypeCat);
% [sigPairs, statsTbl, lme, workTbl, dodTbl] = runLMMStat(T_Vin, T_Vin.FA_Count, ...
%     'CombinedVec',             cellstr(combo), ...
%     'BioRepCol',               'expCat', ...
%     'ContrastExcludePatterns', Input.ContrastExcludePatterns, ...
%     'ContrastRemove',          Input.ContrastRemove, ...
%     'DiffOfDiffs',             'auto');
combo = string(T_VECad.diseaseCat) + "_" + string(T_VECad.extrainfoCat) + "_" + string(T_VECad.sampletypeCat);
[sigPairs, statsTbl, lme, workTbl, dodTbl] = runLMMStat(T_VECad, T_VECad.Av_AverageIntensity, ...
    'CombinedVec',             cellstr(combo), ...
    'BioRepCol',               'expCat', ...
    'ContrastExcludePatterns', Input.ContrastExcludePatterns, ...
    'ContrastRemove',          Input.ContrastRemove, ...
    'DiffOfDiffs',             'auto');
% Genotype x stiffness (Fig 4)
% [intTbl_kPa, ctrTbl_kPa] = runLMMInteractionStat(T_Vin, T_Vin.FA_Count, ...
%     'FactorLevels', {'1kPa','10kPa'});
% % Genotype x ROCK (Fig 3)
% [intTbl_ROCK, ctrTbl_ROCK] = runLMMInteractionStat(T_Vin, T_Vin.FA_Count, ...
%     'FactorLevels', {'rest','R1','R2'});

% ============================================================
%% ============ 6_IMAGES =================
% =============================================================
%saveVectorFigure(gcf, 'C:\Users\sarah\OneDrive - KU Leuven\Analysis_Codes_Matlab\Code saveguarding Jan26\20260114_FDA_V3_SarahBulk\0_PAPER IMAGES\
% VinCTCCM_Gel', 'svg')

%% Overview chronologically
% h = PlotFRETTSCCMBoxplot(T_Vin, 26, 'FRETav', {}, {'R1', 'R2', 'kPa'}); 
% h = PlotFRETTSCCMBoxplot(T_Vin, 26, 'FRETav', {}, {'R1', 'R2'}); 
% h = PlotFRETTSCCMBoxplot(T_Vin, 26, 'FRETav', {}, {'kPa'});
% [T_Vin, BaselineSummary] = Norm_ToCTBaseline(T_Vin, 'FRETav', 1, false);
% h = PlotFRETTSCCMBoxplot(T_Vin, 26, 'FRETavNorm', {'VinTS'}, {'R1', 'R2', 'kPa', 'si48'});
% h = PlotFRETTSCCMBoxplot(T_VECad, 31, 'FRETav', {}, {'R1', 'R2', 'kPa'});
% [T_VECad, BaselineSummary] = Norm_ToCTBaseline(T_VECad, 'FRETav', 1, false);
% h = PlotFRETTSCCMBoxplot(T_VECad, 31, 'FRETavNorm', {'VECadTS'}, {'R1', 'R2', 'kPa', 'si48'});


%% !! Additional
Input.includePatterns_VinFL = {'VinTS_FL', 'VinTL_FL'}; % CAREFUL! : include with ... ! or ! with ...
Input.excludePatterns_Vin = {'ish', 'si48', '20250922','20251116'};
T_Vin = Filter_IncludeExclude(T_Vin, Input.includePatterns_VinFL, Input.excludePatterns_Vin);
%%
Input.includePatterns_VECadFL = {'VECadTS_FL', 'VECadTL_FL'};
Input.excludePatterns_VECad = {'ish', 'si48', '20250922', '20251116'};
T_VECad = Filter_IncludeExclude(T_VECad, Input.includePatterns_VECadFL, Input.excludePatterns_VECad);

% ====================================================================
%% ----------------- FIG_2: Vin on glass  -----------------
% ====================================================================
%% 6 Pannel
PaperImg_ExportIntFRETimage(T_Vin, Input, '20260418_VinTS_FL_CT_si72', 9, ...
    FRET=false, Int=false, Class=true, IntLow = 0, IntHigh = 99.99, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_Vin, Input, '20260418_VinTS_FL_CCM_si72', 12, ...
    FRET=false, Int=false, Class=true, IntLow = 0, IntHigh = 99.5, ExportSVG = true)
%% Intensity, Count, FRET
[~, ~, overviewTbl] = plotMetricScatter(T_Vin, 'FA_Count', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ... % '20250922'
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%display(overviewTbl);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_VinCTCCM_Glass_Count'), 'svg', Width=10, Height=13);
plotMetricScatter(T_Vin, 'Av_AlignIndex', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_VinCTCCM_Glass_Align'), 'svg', Width=10, Height=13);
plotMetricScatter(T_Vin, 'Av_AverageIntensity', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_VinCTCCM_Glass_AvInt'), 'svg', Width=10, Height=13);
plotMetricScatter(T_Vin, 'FRETav', levels={'rest'}, excludePatterns={}, YLimMinSet=15, YLimMaxSet=35, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_VinCTCCM_Glass_FRETav'), 'svg', Width=10, Height=13);
%% Stacked Histogram
Plot_ShapeClassStackedHistogram(T_Vin, Levels={'rest'}, excludePatterns={'TL'}, legendNames = Input.CatNames)
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_VinCTCCM_Glass_ClassHist'), 'svg');

%% SI:
%% Figure S1: Spring choice (T_VinSpring) - ttest!!
% ! NEED TO ADJUST in Illustrator: delete meanBioRep, change color scatter
plotMetricScatter(T_VinSpring, 'FRETav', excludePatterns={'CCM'}, YLimMinSet=15, XCat='springCat', compareCat='sampletypeCat', UseThirdCat=false,...
    Contrasts={{'TS_F40', 'TL_F40'},{'TS_FL', 'TL_FL'},{'TS_stHP35', 'TL_stHP35'}}, CollapseRep=false,...
    Unit='(%) - Apparent', ConditionOrder=Input.SpringOrder, statTest2='ttest', TTestType='independent',...
    paletteName='blues');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTSpring_Glass'), 'svg', Width=12, Height=14);
%% Figure S3: Fraction, intensity, FRET per classification
%% 
plotMetricScatter(T_Vin, T_Vin.FA_summaryFrac(:,1), VarNaming = 'Linear', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summaryFrac', VarColIndex=1, YLimMaxSet = 0.6);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_Glass_Class1Frac'), 'svg', Width=8, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFrac(:,2), VarNaming = 'Elliptical', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summaryFrac', VarColIndex=2, YLimMaxSet = 0.6);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_Glass_Class2Frac'), 'svg', Width=8, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFrac(:,3), VarNaming = 'Small Circular', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='yellows', VarColName='FA_summaryFrac', VarColIndex=3, YLimMaxSet = 0.6);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_Glass_Class3Frac'), 'svg', Width=8, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFrac(:,4), VarNaming = 'Large Circular', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='purples', VarColName='FA_summaryFrac', VarColIndex=4, YLimMaxSet = 0.6);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_Glass_Class4Frac'), 'svg', Width=8, Height=12);
%% 
plotMetricScatter(T_Vin, T_Vin.FA_summaryIntPx(:,1), VarNaming = 'Linear', levels={'rest'}, excludePatterns={}, YLimMinSet=0, YLimMaxSet = 70,...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summaryIntPx', VarColIndex=1);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_Glass_Class1Int'), 'svg', Width=8, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryIntPx(:,2), VarNaming = 'Elliptical', levels={'rest'}, excludePatterns={}, YLimMinSet=0, YLimMaxSet = 70,...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summaryIntPx', VarColIndex=2);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_Glass_Class2Int'), 'svg', Width=8, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryIntPx(:,3), VarNaming = 'Small Circular', levels={'rest'}, excludePatterns={}, YLimMinSet=0, YLimMaxSet = 70,...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='yellows', VarColName='FA_summaryIntPx', VarColIndex=3);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_Glass_Class3Int'), 'svg', Width=8, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryIntPx(:,4), VarNaming = 'Large Circular', levels={'rest'}, excludePatterns={}, YLimMinSet=0, YLimMaxSet = 70,...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='purples', VarColName='FA_summaryIntPx', VarColIndex=4);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_Glass_Class4Int'), 'svg', Width=8, Height=12);
%% 
plotMetricScatter(T_Vin, T_Vin.FA_summaryFRETPx(:,1), VarNaming = 'Linear', levels={'rest'}, excludePatterns={}, YLimMinSet=15, YLimMaxSet = 35,...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summaryFRETPx', VarColIndex=1);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_Glass_Class1FRET'), 'svg', Width=8, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFRETPx(:,2), VarNaming = 'Elliptical', levels={'rest'}, excludePatterns={}, YLimMinSet=15, YLimMaxSet = 35,...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summaryFRETPx', VarColIndex=2);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_Glass_Class2FRET'), 'svg', Width=8, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFRETPx(:,3), VarNaming = 'Small Circular', levels={'rest'}, excludePatterns={}, YLimMinSet=15, YLimMaxSet = 35,...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='yellows', VarColName='FA_summaryFRETPx', VarColIndex=3);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_Glass_Class3FRET'), 'svg', Width=8, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFRETPx(:,4), VarNaming = 'Large Circular', levels={'rest'}, excludePatterns={}, YLimMinSet=15, YLimMaxSet = 35,...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='purples', VarColName='FA_summaryFRETPx', VarColIndex=4);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_Glass_Class4FRET'), 'svg', Width=8, Height=12);
%% Figure S4: Fraction, intensity, FRET between classification
% Comparing categories
longTblFRET = Help_expandMatrixColumn(T_Vin, 'FA_summaryFRETPx', Input.CatNames);
longTblInt = Help_expandMatrixColumn(T_Vin, 'FA_summaryIntPx', Input.CatNames);
%%
plotMetricScatter(longTblInt, 'FA_summaryIntPx', excludePatterns={'VinTL', 'CCM', 'R1', 'R2', 'kPa'}, YLimMinSet=0, YLimMaxSet = 90,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='rest',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCT_Intcat'), 'svg', Width=13, Height=13);
plotMetricScatter(longTblInt, 'FA_summaryIntPx', excludePatterns={'VinTL', 'CT', 'R1', 'R2', 'kPa'}, YLimMinSet=0, YLimMaxSet = 90,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='rest',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCCM_Intcat'), 'svg', Width=13, Height=13);
%%
plotMetricScatter(longTblFRET, 'FA_summaryFRETPx', excludePatterns={'VinTL', 'CCM', 'R1', 'R2', 'kPa'}, YLimMinSet=15, YLimMaxSet = 38,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='rest',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCT_FRETcat'), 'svg', Width=13, Height=13);
plotMetricScatter(longTblFRET, 'FA_summaryFRETPx', excludePatterns={'VinTL', 'CT', 'R1', 'R2', 'kPa'}, YLimMinSet=15, YLimMaxSet = 38,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='rest',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCCM_FRETcat'), 'svg', Width=13, Height=13);
%% Figure S5: TL CTvsCCM
plotMetricScatter(T_Vin, 'FRETav', levels={'rest'}, excludePatterns={'TS'}, YLimMinSet=15, YLimMaxSet=38,...
    CondColors=Input.condColors, statTest2='ttest', TTestType='independent', ThirdCatFilter={'TL'}, Contrasts={{'CT_rest_TL', 'CCM_rest_TL'}});
saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_TL_Glass_FRETav'), 'svg', Width=10, Height=13);
%% SI - extra
% %% Plot morph vs FRET / FRET vs Int
% Plot_CompareVar(T_Vin, "AverageIntensity", "FRET", ...
%     Levels       = {'rest', 'R1', 'R2', '10kPa', '1kPa'}, ...
%     PlotType     = 'scatter', ...
%     ScatterMode  = 'roi',...       % ScatterMode = 'all',...
%     ScatterStyle = 'points',...
%     CondColors   = Input.condColors);
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_Vin_AllStat_vsFRET_scatter'), 'svg');
% %%
% Plot_CompareVar(T_Vin, "AverageIntensity", "FRET", ...
%     Levels = {'rest', 'R1', 'R2', '10kPa', '1kPa'}, PlotType = 'binned', CondColors = Input.condColors);
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_Vin_AllStat_vsFRET_ribbon'), 'svg');
% %% FL choice - TL CTvsCCM - more
% % % % Optional: add morphology analysis for TLs as well?
% % plotMetricScatter(T_Vin, 'FA_Count', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
% %     Unit='Fraction', ContrastExcludePatterns = Input.ContrastExcludePatterns,...
% %     CondColors=Input.condColors, ThirdCatFilter={'TL'});
% % %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_TL_Glass_FACount'), 'svg', Width=8, Height=14);
% % plotMetricScatter(T_Vin, 'Av_AverageIntensity', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
% %     Unit='', ContrastExcludePatterns = Input.ContrastExcludePatterns,...
% %     CondColors=Input.condColors, ThirdCatFilter={'TL'});
% % %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_TL_Glass_AvInt'), 'svg', Width=8, Height=14);
% % plotMetricScatter(T_Vin, T_Vin.FA_summaryFrac(:,1), VarNaming = 'Linear', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
% %     Unit='Fraction', ContrastExcludePatterns = Input.ContrastExcludePatterns,...
% %     paletteName='bluepurples', VarColName='FA_summaryFrac', VarColIndex=1, YLimMaxSet = 0.6, ThirdCatFilter={'TL'});
% % %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_TL_Glass_Class1'), 'svg', Width=8, Height=14);
% % plotMetricScatter(T_Vin, T_Vin.FA_summaryFrac(:,4), VarNaming = 'Large Circular', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
% %     Unit='Fraction', ContrastExcludePatterns = Input.ContrastExcludePatterns,...
% %     paletteName='purples', VarColName='FA_summaryFrac', VarColIndex=4, YLimMaxSet = 0.6, ThirdCatFilter={'TL'});
% % %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_VinCTCCM_TL_Glass_Class4'), 'svg', Width=8, Height=14);
% %% Statistics: where do the stats come from (2 orderings):
% plotMetricScatter(T_Vin, 'FRETav', levels={'rest', 'R1', 'R2', '10kPa', '1kPa'}, excludePatterns={}, YLimMinSet=15,...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, ConditionOrder=Input.CondOrder,...
%     CondColors=Input.condColors, ThirdCatFilter={'TS', 'TL'});
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_Vin_AllStat_Order1'), 'svg', Width=12, Height=14);
% plotMetricScatter(T_Vin, 'FRETav', levels={'rest', 'R1', 'R2', '10kPa', '1kPa'}, excludePatterns={}, YLimMinSet=15,...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, ConditionOrder=Input.CondOrder2,...
%     CondColors=Input.condColors, ThirdCatFilter={'TS', 'TL'});
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_Vin_AllStat_Order2'), 'svg', Width=12, Height=14);
% %% Live vs fixed
% plotMetricScatter(T_VinLiveFixed, 'FRETav', levels={'rest'}, excludePatterns={}, YLimMinSet=15, CollapseRep=false, ThirdCat='LiveFixedCat',...
%     ThirdCatFilter={'Live','Fixed'}, CondColors=Input.condColors);
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig2_SI_Vin_LiveFixed'), 'svg', Width=8, Height=14);
% %% Other code
% % Plot_DistrubutionVar_ByCond(T_VECad, "CellSize", 10);
% % Plot_FAIntClassComparison(T_Vin, XLimit=400, excludePatterns={'TL'});
% % Plot_FAClassDistribution(T_Vin, Metric='FRET', Levels={'rest', '10kPa', '1kPa'}, excludePatterns={'si48', 'ish', 'TL'}); % FRET/Int/Frac
% % Plot_ShapeClassEvolution(T_Vin, 'FRET', excludePatterns={'si48', 'ish', 'TL'}, Levels={'rest', '10kPa', '1kPa'}); % FRET/Int/Frac
% % Plot_ColorFAMask(T_Vin, '20251117_VinTS_FL_CCM_si72_24h', 4)
% % Plot_VisualizeOrientation(T_Vin, 127);
% % plotFRETchronological(T_Vin); % check if chronological has effect
% % mask = T_Vin.FRETindex{i}; QC = T_Vin.('Area'){i}; QC_FA_Diagnostic(mask, QC);
% % plotScatter2var(T_Vin, 'FRETav', 'Av_AverageIntensity', 'excludePatterns', {'TL'}, 'XCutoff', 15, 'YCutoff', 0, 'Alpha', 0.5);
% 
% %Plot_DistrubutionVar_ByCond(T_Vin, "Area", nBins=100, filterVar="MinMajRatio", filterOp=">=", filterVal=0.5);
% %Plot_DistrubutionVar_ByCond(T_Vin, "MinMajRatio", nBins=100);%, filterVar="MinMajRatio", filterOp="<", filterVal=0.5);
% %Plot_DistrubutionVar_ByCond(T_Vin, "MajAx", nBins=100);
% Plot_DistrubutionVar_ByCond(T_Vin, "MajAx", nBins=100, filterVar1="MinMajRatio", filterOp1=">=", filterVal1=0.5,...
%     filterVar2="Area", filterOp2=">", filterVal2=25);
% %Plot_DistrubutionVar_ByCond(T_Vin, "MajAx", nBins=100, filterVar1="MinMajRatio", filterOp1="<", filterVal1=0.5);


% ====================================================================

%% ----------------- FIG_3: VECad on glass  -----------------
% ====================================================================
%% 6 Pannel
PaperImg_ExportIntFRETimage(T_VECad, Input, '20260414_VECadTS_FL_CT_si72', 8, ...
    Class=true, Skel=true, ThickMax=11, IntLow = 0, IntHigh = 99.99, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_VECad, Input, '20251118_VECadTS_FL_CCM_si72', 2, ...
    Class=true, Skel=true, ThickMax=11, IntLow = 0, IntHigh = 99.5, ExportSVG = true)
%% CellShape, Thickness(skel), Intensity, FRET
plotMetricScatter(T_VECad, 'CellAv_MinMajRatio', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_VECadCTCCM_Glass_CellShape'), 'svg', Width=10, Height=13);
plotMetricScatter(T_VECad, 'Thickness_Skeleton_Av', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_VECadCTCCM_Glass_Thick'), 'svg', Width=10, Height=13);
plotMetricScatter(T_VECad, 'Av_AverageIntensity', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
    Unit='', ContrastExcludePatterns = Input.ContrastExcludePatterns, ...
    CondColors=Input.condColors, ThirdCatFilter={'TS'});
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_VECadCTCCM_Glass_AvInt'), 'svg', Width=10, Height=13);
plotMetricScatter(T_VECad, 'FRETav', levels={'rest'}, excludePatterns={}, YLimMinSet=15, Unit='(%) - Apparent',...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors, YLimMaxSet=40);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_VECadCTCCM_Glass_FRETav'), 'svg', Width=10, Height=13);
%% Histogram
Plot_ShapeClassStackedHistogram(T_VECad, Levels={'rest'}, excludePatterns={'TL'}, legendNames = Input.CatNamesAJ)
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_VECadCTCCM_Glass_ClassHist'), 'svg');

%% SI
%% Figure S1: Spring choice (T_VECadSpring) - ttest!!
plotMetricScatter(T_VECadSpring, 'FRETav', excludePatterns={'CCM'}, YLimMinSet=15, XCat='springCat', compareCat='sampletypeCat', UseThirdCat=false,...
    Contrasts={{'TS_F40', 'TL_F40'},{'TS_FL', 'TL_FL'},{'TS_stHP35', 'TL_stHP35'}}, CollapseRep=false,...
    Unit='(%) - Apparent', ConditionOrder=Input.SpringOrder, statTest2='ttest', TTestType='independent',...
    paletteName='blues');
saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCTSpring_Glass'), 'svg', Width=12, Height=14);
%% Figure S7: Cell size
plotMetricScatter(T_VECad, 'CellAvUmSize', levels={'rest'}, excludePatterns={}, YLimMinSet=0, Unit='(µm)', ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCTCCM_Glass_CellSize'), 'svg', Width=8, Height=12);
%% Figure S9: Fraction, intensity, FRET per classification
plotMetricScatter(T_VECad, T_VECad.FA_summaryFrac(:,1), VarNaming = 'Reticular', levels={'rest'}, excludePatterns={}, YLimMinSet=0, YLimMaxSet = 1.19,...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summaryFrac', VarColIndex=1);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCTCCM_Glass_Class1Frac'), 'svg', Width=8, Height=12);
plotMetricScatter(T_VECad, T_VECad.FA_summaryFrac(:,2), VarNaming = 'Linear', levels={'rest'}, excludePatterns={}, YLimMinSet=0, YLimMaxSet = 1.19,...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summaryFrac', VarColIndex=2);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCTCCM_Glass_Class2Frac'), 'svg', Width=8, Height=12);
%%
plotMetricScatter(T_VECad, T_VECad.FA_summarySkelInt(:,1), VarNaming = 'Reticular', levels={'rest'}, excludePatterns={}, YLimMinSet=0, YLimMaxSet = 90,...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summarySkelInt', VarColIndex=1);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCTCCM_Glass_Class1Int'), 'svg', Width=8, Height=12);
plotMetricScatter(T_VECad, T_VECad.FA_summarySkelInt(:,2), VarNaming = 'Linear', levels={'rest'}, excludePatterns={}, YLimMinSet=0, YLimMaxSet = 90,...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summarySkelInt', VarColIndex=2);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCTCCM_Glass_Class2Int'), 'svg', Width=8, Height=12);
%%
plotMetricScatter(T_VECad, T_VECad.FA_summarySkelFRET(:,1), VarNaming = 'Reticular', levels={'rest'}, excludePatterns={}, YLimMinSet=15, YLimMaxSet = 40,...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summarySkelFRET', VarColIndex=1);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCTCCM_Glass_Class1FRET'), 'svg', Width=8, Height=12);
plotMetricScatter(T_VECad, T_VECad.FA_summarySkelFRET(:,2), VarNaming = 'Linear', levels={'rest'}, excludePatterns={}, YLimMinSet=15, YLimMaxSet = 40,...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summarySkelFRET', VarColIndex=2);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCTCCM_Glass_Class2FRET'), 'svg', Width=8, Height=12);
%% Figure S10: Fraction, intensity, FRET between classification
longTblFRETVE = Help_expandMatrixColumn(T_VECad, 'FA_summarySkelFRET', Input.CatNamesAJ);
longTblIntVE = Help_expandMatrixColumn(T_VECad, 'FA_summarySkelInt', Input.CatNamesAJ);
longTblThickVE = Help_expandMatrixColumn(T_VECad, 'FA_summarySkelThickness', Input.CatNamesAJ);
%%
plotMetricScatter(longTblIntVE, 'FA_summarySkelInt', excludePatterns={'VECadTL', 'CCM', 'R1', 'R2', 'kPa'}, YLimMinSet=0, YLimMaxSet = 90,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='rest',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCT_Intcat'), 'svg', Width=13, Height=13);
plotMetricScatter(longTblIntVE, 'FA_summarySkelInt', excludePatterns={'VECadTL', 'CT', 'R1', 'R2', 'kPa'}, YLimMinSet=0, YLimMaxSet = 90,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='rest',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCCM_Intcat'), 'svg', Width=13, Height=13);
%%
plotMetricScatter(longTblFRETVE, 'FA_summarySkelFRET', excludePatterns={'VECadTL', 'CCM', 'R1', 'R2', 'kPa'}, YLimMinSet=15, YLimMaxSet = 40,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='rest',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCT_FRETcat'), 'svg', Width=13, Height=13);
plotMetricScatter(longTblFRETVE, 'FA_summarySkelFRET', excludePatterns={'VECadTL', 'CT', 'R1', 'R2', 'kPa'}, YLimMinSet=15, YLimMaxSet = 40,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='rest',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCCM_FRETcat'), 'svg', Width=13, Height=13);
%%  Figure S11: TL CTvsCCM
plotMetricScatter(T_VECad, 'FRETav', levels={'rest'}, excludePatterns={'TS'}, YLimMinSet=15, YLimMaxSet=38,...
    CondColors=Input.condColors, statTest2='ttest', TTestType='independent', ThirdCatFilter={'TL'}, Contrasts={{'CT_rest_TL', 'CCM_rest_TL'}});
saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCTCCM_TL_Glass_FRETav'), 'svg', Width=10, Height=13);
%% SI - extra
% %% Plot morph vs FRET / FRET vs Int STILL TO DO
% Plot_CompareVar(T_VECad, "Thickness_Skeleton", "Skel_FRET", ...
%     Levels       = {'rest', 'R1', 'R2', '10kPa', '1kPa'}, ...
%     PlotType     = 'scatter', ...
%     ScatterMode  = 'roi',...       % ScatterMode = 'all',...
%     ScatterStyle = 'points',...
%     CondColors   = Input.condColors);
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECad_AllStat_ThickvsFRET_scatter'), 'svg');
% %%
% Plot_CompareVar(T_VECad, "Thickness_Skeleton", "Skel_FRET", ...
%     Levels = {'rest', 'R1', 'R2', '10kPa', '1kPa'}, PlotType = 'binned', CondColors = Input.condColors);
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECad_AllStat_ThickvsFRET_ribbon'), 'svg');
% %%
% Plot_CompareVar(T_VECad, "Thickness_Skeleton", "Skel_Int", ...
%     Levels = {'rest', 'R1', 'R2', '10kPa', '1kPa'}, PlotType = 'binned', CondColors = Input.condColors);
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECad_AllStat_ThickvsInt_ribbon'), 'svg');
% Plot_CompareVar(T_VECad, "Skel_Int", "Skel_FRET", ...
%     Levels = {'rest', 'R1', 'R2', '10kPa', '1kPa'}, PlotType = 'binned', CondColors = Input.condColors);
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECad_AllStat_IntvsFRET_ribbon'), 'svg');
% %% FL choice - TL CTvsCCM
% % ! NEED TO ADJUST in Illustrator: delete meanBioRep, change color scatter
% % paired??
% plotMetricScatter(T_VECad, 'FRETav', levels={'rest'}, excludePatterns={}, YLimMinSet=15, ...
%     Unit='(%) - Apparent', ContrastExcludePatterns = Input.ContrastExcludePatterns, ...
%     CondColors=Input.condColors, ThirdCatFilter={'TL'});
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCTCCM_TL_Glass_FRETav'), 'svg', Width=8, Height=14);
% % % % Optional: add morphology analysis for TLs as well?
% % plotMetricScatter(T_VECad, 'CellAvUmSize', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
% %     Unit='(µm)', ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors, ThirdCatFilter={'TL'});
% % %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCTCCM_TL_Glass_CellSize'), 'svg', Width=8, Height=14);
% % plotMetricScatter(T_VECad, 'CellAv_MinMajRatio', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
% %     Unit='', ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors, ThirdCatFilter={'TL'});
% % %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCTCCM_TL_Glass_CellShape'), 'svg', Width=8, Height=14);
% % plotMetricScatter(T_VECad, 'Thickness_Skeleton_Av', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
% %     Unit='(µm?)', ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors, ThirdCatFilter={'TL'});
% % saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCTCCM_TL_Glass_Thick'), 'svg', Width=8, Height=14);
% % plotMetricScatter(T_VECad, 'Av_AverageIntensity', levels={'rest'}, excludePatterns={}, YLimMinSet=0, ...
% %     Unit='', ContrastExcludePatterns = Input.ContrastExcludePatterns, ...
% %     CondColors=Input.condColors, ThirdCatFilter={'TL'});
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECadCTCCM_TL_Glass_AvInt'), 'svg', Width=8, Height=14);
% %% Statistics: where do the stats come from (2 orderings):
% plotMetricScatter(T_VECad, 'FRETav', levels={'rest', 'R1', 'R2', '10kPa', '1kPa'}, excludePatterns={}, YLimMinSet=15,...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, ConditionOrder=Input.CondOrder,...
%     CondColors=Input.condColors, ThirdCatFilter={'TS', 'TL'});
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECad_AllStat_Order1'), 'svg', Width=16, Height=14);
% plotMetricScatter(T_VECad, 'FRETav', levels={'rest', 'R1', 'R2', '10kPa', '1kPa'}, excludePatterns={}, YLimMinSet=15,...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, ConditionOrder=Input.CondOrder2,...
%     CondColors=Input.condColors, ThirdCatFilter={'TS', 'TL'});
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECad_AllStat_Order2'), 'svg', Width=16, Height=14);
% %% Live vs fixed
% plotMetricScatter(T_VECadLiveFixed, 'FRETav', levels={'rest'}, excludePatterns={}, YLimMinSet=15, CollapseRep=false, ThirdCat='LiveFixedCat',...
%     ThirdCatFilter={'Live','Fixed'}, CondColors=Input.condColors);
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig3_SI_VECad_LiveFixed'), 'svg', Width=8, Height=14);
% %% Other code
% % %% Distribution - Thickness - Skeleton
% % [G, fig] = plotThickness(T_VECad, 'Thickness_Skeleton', 0.5, 10, {'VECadTS_FL'}, {'24h','si48','ish', 'TL', '202508','202509'});
% % %% DistributionCumm - Thickness - Skeleton
% % [G, fig] = plotThicknessCumulative(T_VECad, 'Thickness_Skeleton', Input, excludePatterns={'si48', 'ish','TL','202508','202509','20251124'});
% % Plot_ColorCellCellMask(T_VECad, '20251118_VECadTS_FL_CCM_si72', 2)
% % plotScatter2var(T_VECad, 'Thickness_Skeleton_Av', 'FRETav', 'excludePatterns', {'TL'}, 'XCutoff', 0, 'YCutoff', 0, 'Alpha', 0.5);
% % plotFRETvsThickness(T_VECad, 'XVar','thickness', 'YVar','fret','excludePatterns', {'si48', 'ish','TL','20251124'}); % fret/int/thickness


% ====================================================================
%% ----------------- FIG_4: Vin with ROCK  -----------------
% ====================================================================
%% 18 Pannel (6 sets)
PaperImg_ExportIntFRETimage(T_Vin, Input, '20260427_VinTS_FL_CT_glass_si72', 16, ...
    Class=true, IntLow = 0, IntHigh = 99.99, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_Vin, Input, '20260427_VinTS_FL_CCM_glass_si72', 7, ...
    Class=true, IntLow = 0, IntHigh = 99.5, ExportSVG = true)

PaperImg_ExportIntFRETimage(T_Vin, Input, '20260415_VinTS_FL_CT_R1_si72', 9, ...
    Class=true, IntLow = 0, IntHigh = 99.99, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_Vin, Input, '20260415_VinTS_FL_CT_R2_si72', 14, ...
    Class=true, IntLow = 0, IntHigh = 99.5, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_Vin, Input, '20260415_VinTS_FL_CCM_R1_si72', 4, ...
    Class=true, IntLow = 0, IntHigh = 99.99, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_Vin, Input, '20260418_VinTS_FL_CCM_R2_si72', 11, ...
    Class=true, IntLow = 0, IntHigh = 99.5, ExportSVG = true)
%% Intensity, Count, FRET - no CT_R1 and CT_R2
plotMetricScatter(T_Vin, 'FA_Count', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0,  ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_VinCTCCM_ROCK_Count_filt'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, 'Av_AlignIndex', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_VinCTCCM_Glass_Align_filt'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, 'Av_AverageIntensity', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0,  ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_VinCTCCM_ROCK_AvInt_filt'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, 'FRETav', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=15,  YLimMaxSet=35, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_VinCTCCM_ROCK_FRETav_filt'), 'svg', Width=16, Height=12);
%% Histogram
Plot_ShapeClassStackedHistogram(T_Vin, Levels={'rest', 'R1', 'R2'}, excludePatterns={'TL'}, legendNames = Input.CatNames)
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_VinCTCCM_ROCK_ClassHist'), 'svg');

%% SI:
%% Figure S12: Main figure with all ROCKs
plotMetricScatter(T_Vin, 'FA_Count', levels={'rest', 'R1', 'R2'}, excludePatterns={}, YLimMinSet=0,  ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_Count'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, 'Av_AlignIndex', levels={'rest', 'R1', 'R2'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_Glass_Align'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, 'Av_AverageIntensity', levels={'rest', 'R1', 'R2'}, excludePatterns={}, YLimMinSet=0,  ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_AvInt'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, 'FRETav', levels={'rest', 'R1', 'R2'}, excludePatterns={}, YLimMinSet=15, YLimMaxSet=35, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_FRETav'), 'svg', Width=16, Height=12);
%% Figure S13: Fraction, intensity, FRET per classification
plotMetricScatter(T_Vin, T_Vin.FA_summaryFrac(:,1), VarNaming = 'Linear', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summaryFrac', VarColIndex=1, YLimMaxSet = 0.7);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_Class1'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFrac(:,2), VarNaming = 'Elliptical', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summaryFrac', VarColIndex=2, YLimMaxSet = 0.7);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_Class2'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFrac(:,3), VarNaming = 'Small Circular', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='yellows', VarColName='FA_summaryFrac', VarColIndex=3, YLimMaxSet = 0.7);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_Class3'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFrac(:,4), VarNaming = 'Large Circular', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='purples', VarColName='FA_summaryFrac', VarColIndex=4, YLimMaxSet = 0.7);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_Class4'), 'svg', Width=16, Height=12);
%% 
plotMetricScatter(T_Vin, T_Vin.FA_summaryIntPx(:,1), VarNaming = 'Linear', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summaryIntPx', VarColIndex=1, YLimMaxSet = 75);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_Class1Int'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryIntPx(:,2), VarNaming = 'Elliptical', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summaryIntPx', VarColIndex=2, YLimMaxSet = 75);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_Class2Int'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryIntPx(:,3), VarNaming = 'Small Circular', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='yellows', VarColName='FA_summaryIntPx', VarColIndex=3, YLimMaxSet = 75);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_Class3Int'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryIntPx(:,4), VarNaming = 'Large Circular', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='purples', VarColName='FA_summaryIntPx', VarColIndex=4, YLimMaxSet = 75);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_Class4Int'), 'svg', Width=16, Height=12);
%% 
plotMetricScatter(T_Vin, T_Vin.FA_summaryFRETPx(:,1), VarNaming = 'Linear', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=15, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summaryFRETPx', VarColIndex=1, YLimMaxSet = 35);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_Class1FRET'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFRETPx(:,2), VarNaming = 'Elliptical', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=15, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summaryFRETPx', VarColIndex=2, YLimMaxSet = 35);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_Class2FRET'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFRETPx(:,3), VarNaming = 'Small Circular', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=15, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='yellows', VarColName='FA_summaryFRETPx', VarColIndex=3, YLimMaxSet = 35);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_Class3FRET'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFRETPx(:,4), VarNaming = 'Large Circular', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=15, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='purples', VarColName='FA_summaryFRETPx', VarColIndex=4, YLimMaxSet = 35);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCTCCM_ROCK_Class4FRET'), 'svg', Width=16, Height=12);
%% Figure S14: Intensity, FRET between classification
% plotMetricScatter(longTblInt, 'FA_summaryIntPx', excludePatterns={'CCM'}, includePatterns={'R1'}, YLimMinSet=0,...
%     XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R1',...
%     CondColors=Input.condColorsClass, statTest2='rm');
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCT_R1_Intcat'), 'svg', Width=16, Height=12);
% plotMetricScatter(longTblInt, 'FA_summaryIntPx', excludePatterns={'CCM'}, includePatterns={'R2'}, YLimMinSet=0,...
%     XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R2',...
%     CondColors=Input.condColorsClass, statTest2='rm');
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCT_R2_Intcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblInt, 'FA_summaryIntPx', excludePatterns={'CT'}, includePatterns={'R1'}, YLimMinSet=0, YLimMaxSet = 75,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R1',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCCM_R1_Intcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblInt, 'FA_summaryIntPx', excludePatterns={'CT'}, includePatterns={'R2'}, YLimMinSet=0, YLimMaxSet = 75,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R2',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCCM_R2_Intcat'), 'svg', Width=16, Height=12);
%%
% plotMetricScatter(longTblFRET, 'FA_summaryFRETPx', excludePatterns={'CCM'}, includePatterns={'R1'}, YLimMinSet=15,...
%     XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R1',...
%     CondColors=Input.condColorsClass, statTest2='rm');
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCT_R1_FRETcat'), 'svg', Width=16, Height=12);
% plotMetricScatter(longTblFRET, 'FA_summaryFRETPx', excludePatterns={'CCM'}, includePatterns={'R2'}, YLimMinSet=15,...
%     XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R2',...
%     CondColors=Input.condColorsClass, statTest2='rm');
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCT_R2_FRETcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblFRET, 'FA_summaryFRETPx', excludePatterns={'CT'}, includePatterns={'R1'}, YLimMinSet=15, YLimMaxSet = 36,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R1',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCCM_R1_FRETcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblFRET, 'FA_summaryFRETPx', excludePatterns={'CT'}, includePatterns={'R2'}, YLimMinSet=15, YLimMaxSet = 36,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R2',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig4_SI_VinCCM_R2_FRETcat'), 'svg', Width=16, Height=12);


% ====================================================================
%% ----------------- FIG_5: VECad with ROCK  -----------------
% ====================================================================
%% 18 Pannel (6 sets)
PaperImg_ExportIntFRETimage(T_VECad, Input, '20260427_VECadTS_FL_CT_glass_si72', 5, ...
    Class=true, Skel=true, ThickMax=13, IntLow = 0, IntHigh = 99.99, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_VECad, Input, '20260420_VECadTS_FL_CCM_si72', 16, ...
    Class=true, Skel=true, ThickMax=13, IntLow = 0, IntHigh = 99.5, ExportSVG = true)

PaperImg_ExportIntFRETimage(T_VECad, Input, '20260414_VECadTS_FL_CT_R1_si72', 14, ...
    Class=true, Skel=true, ThickMax=13, IntLow = 0, IntHigh = 99.99, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_VECad, Input, '20260414_VECadTS_FL_CT_R2_si72', 11, ...
    Class=true, Skel=true, ThickMax=13, IntLow = 0, IntHigh = 99.5, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_VECad, Input, '20260420_VECadTS_FL_CCM_R1_si72', 7, ...
    Class=true, Skel=true, ThickMax=13, IntLow = 0, IntHigh = 99.99, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_VECad, Input, '20260420_VECadTS_FL_CCM_R2_si72', 11, ...
    Class=true, Skel=true, ThickMax=13, IntLow = 0, IntHigh = 99.5, ExportSVG = true)
%% CellSize, CellShape, Thickness(skel), FRET - no CT-R1 or CT-R2
plotMetricScatter(T_VECad, 'CellAv_MinMajRatio', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, Unit='', ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_VECadCTCCM_ROCK_CellShape_filt'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, 'Thickness_Skeleton_Av', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, Unit='(µm)', ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_VECadCTCCM_ROCK_Thick_filt'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, 'Av_AverageIntensity', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_VECadCTCCM_ROCK_AvInt_filt'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, 'FRETav', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=15, Unit='(%) - Apparent', ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_VECadCTCCM_ROCK_FRETav_filt'), 'svg', Width=16, Height=12);
%% Histogram
Plot_ShapeClassStackedHistogram(T_VECad, Levels={'rest', 'R1', 'R2'}, excludePatterns={'TL'}, legendNames = Input.CatNamesAJ)
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_VECadCTCCM_ROCK_ClassHist'), 'svg');

%% SI:
%% Figure S16: Cell size
plotMetricScatter(T_VECad, 'CellAvUmSize', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, Unit='(µm^2)', ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCTCCM_ROCK_CellSize_filt'), 'svg', Width=16, Height=12);
%% Figure S18: Main figure with all ROCKs
plotMetricScatter(T_VECad, 'CellAv_MinMajRatio', levels={'rest', 'R1', 'R2'}, excludePatterns={}, YLimMinSet=0, Unit='', ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCTCCM_ROCK_CellShape'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, 'Thickness_Skeleton_Av', levels={'rest', 'R1', 'R2'}, excludePatterns={}, YLimMinSet=0, Unit='(µm)', ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCTCCM_ROCK_Thick'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, 'Av_AverageIntensity', levels={'rest', 'R1', 'R2'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCTCCM_ROCK_AvInt'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, 'FRETav', levels={'rest', 'R1', 'R2'}, excludePatterns={}, YLimMinSet=15, Unit='(%) - Apparent', ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCTCCM_ROCK_FRETav'), 'svg', Width=16, Height=12);
%% Figure S19: Fraction, intensity, FRET per classification
plotMetricScatter(T_VECad, T_VECad.FA_summaryFrac(:,1), VarNaming = 'Reticular', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, ... %'CT_R1', 'CT_R2'
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summaryFrac', VarColIndex=1, YLimMaxSet = 1.4);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCTCCM_ROCK_Class1Frac'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, T_VECad.FA_summaryFrac(:,2), VarNaming = 'Linear', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summaryFrac', VarColIndex=2, YLimMaxSet = 1.4);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCTCCM_ROCK_Class2Frac'), 'svg', Width=16, Height=12);
%% 
plotMetricScatter(T_VECad, T_VECad.FA_summarySkelInt(:,1), VarNaming = 'Reticular', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summarySkelInt', VarColIndex=1, YLimMaxSet = 85);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCTCCM_ROCK_Class1Int'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, T_VECad.FA_summarySkelInt(:,2), VarNaming = 'Linear', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summarySkelInt', VarColIndex=2, YLimMaxSet = 85);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCTCCM_ROCK_Class2Int'), 'svg', Width=16, Height=12);
%% 
plotMetricScatter(T_VECad, T_VECad.FA_summarySkelFRET(:,1), VarNaming = 'Reticular', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=15, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summarySkelFRET', VarColIndex=1, YLimMaxSet = 49.5);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCTCCM_ROCK_Class1FRET'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, T_VECad.FA_summarySkelFRET(:,2), VarNaming = 'Linear', levels={'rest', 'R1', 'R2'}, excludePatterns={'CT_R1', 'CT_R2'}, YLimMinSet=15, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summarySkelFRET', VarColIndex=2, YLimMaxSet = 49.5);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCTCCM_ROCK_Class2FRET'), 'svg', Width=16, Height=12);
%% Figure S20: Intensity, FRET between classification
% plotMetricScatter(longTblIntVE, 'FA_summarySkelInt', excludePatterns={'CCM'}, includePatterns={'R1'}, YLimMinSet=0,...
%     XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R1',...
%     CondColors=Input.condColorsClassAJ, statTest2='rm');
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCT_R1_Intcat'), 'svg', Width=16, Height=12);
% plotMetricScatter(longTblIntVE, 'FA_summarySkelInt', excludePatterns={'CCM'}, includePatterns={'R2'}, YLimMinSet=0,...
%     XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R2',...
%     CondColors=Input.condColorsClassAJ, statTest2='rm');
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCT_R2_Intcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblIntVE, 'FA_summarySkelInt', excludePatterns={'CT'}, includePatterns={'R1'}, YLimMinSet=0, YLimMaxSet = 80,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R1',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCCM_R1_Intcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblIntVE, 'FA_summarySkelInt', excludePatterns={'CT'}, includePatterns={'R2'}, YLimMinSet=0, YLimMaxSet = 80,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R2',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCCM_R2_Intcat'), 'svg', Width=16, Height=12);
%%
% plotMetricScatter(longTblFRETVE, 'FA_summarySkelFRET', excludePatterns={'CCM'}, includePatterns={'R1'}, YLimMinSet=15,...
%     XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R1',...
%     CondColors=Input.condColorsClassAJ, statTest2='rm');
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCT_R1_FRETcat'), 'svg', Width=16, Height=12);
% plotMetricScatter(longTblFRETVE, 'FA_summarySkelFRET', excludePatterns={'CCM'}, includePatterns={'R2'}, YLimMinSet=15,...
%     XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R2',...
%     CondColors=Input.condColorsClassAJ, statTest2='rm');
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCT_R2_FRETcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblFRETVE, 'FA_summarySkelFRET', excludePatterns={'CT'}, includePatterns={'R1'}, YLimMinSet=15, YLimMaxSet = 38,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R1',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCCM_R1_FRETcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblFRETVE, 'FA_summarySkelFRET', excludePatterns={'CT'}, includePatterns={'R2'}, YLimMinSet=15, YLimMaxSet = 38,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='R2',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig5_SI_VECadCCM_R2_FRETcat'), 'svg', Width=16, Height=12);


% ====================================================================
%% ----------------- FIG_6: Vin on Gel  -----------------
% ====================================================================
%% 12 Pannel (4 sets)
PaperImg_ExportIntFRETimage(T_Vin, Input, '20260429_VinTS_FL_CT_1kPa_si72', 2, ...
    Class=true, IntLow = 0, IntHigh = 99.99, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_Vin, Input, '20260427_VinTS_FL_CT_10kPa_si72', 6, ...
    Class=true, IntLow = 0, IntHigh = 99.5, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_Vin, Input, '20260427_VinTS_FL_CCM_1kPa_si72', 9, ...
    Class=true, IntLow = 0, IntHigh = 99.99, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_Vin, Input, '20260429_VinTS_FL_CCM_10kPa_si72', 15, ...
    Class=true, IntLow = 0, IntHigh = 99.5, ExportSVG = true)
%% Intensity, Count, FRET
plotMetricScatter(T_Vin, 'FA_Count', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0,  ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors,...
    ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_VinCTCCM_Gel_Count'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, 'Av_AlignIndex', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors,...
    ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_VinCTCCM_Gel_Align'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, 'Av_AverageIntensity', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0,  ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors,...
    ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_VinCTCCM_Gel_AvInt'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, 'FRETav', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=15,  ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors,...
    ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_VinCTCCM_Gel_FRETav'), 'svg', Width=16, Height=12);
%% Histogram
Plot_ShapeClassStackedHistogram(T_Vin, Levels={'1kPa', '10kPa'}, excludePatterns={'TL'}, legendNames = Input.CatNames)
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_VinCTCCM_ROCK_ClassHist'), 'svg');

%% SI:
%% Figure S21: Fraction, intensity, FRET per classification
plotMetricScatter(T_Vin, T_Vin.FA_summaryFrac(:,1), VarNaming = 'Linear', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summaryFrac', VarColIndex=1, YLimMaxSet = 0.8, ...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCTCCM_Gel_Class1'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFrac(:,2), VarNaming = 'Elliptical', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summaryFrac', VarColIndex=2, YLimMaxSet = 0.8,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCTCCM_Gel_Class2'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFrac(:,3), VarNaming = 'Small Circular', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='yellows', VarColName='FA_summaryFrac', VarColIndex=3, YLimMaxSet = 0.8,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCTCCM_Gel_Class3'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFrac(:,4), VarNaming = 'Large Circular', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='purples', VarColName='FA_summaryFrac', VarColIndex=4, YLimMaxSet = 0.8,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCTCCM_Gel_Class4'), 'svg', Width=16, Height=12);
%%
plotMetricScatter(T_Vin, T_Vin.FA_summaryIntPx(:,1), VarNaming = 'Linear', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summaryIntPx', VarColIndex=1, YLimMaxSet = 110,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCTCCM_Glass_Class1Int'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryIntPx(:,2), VarNaming = 'Elliptical', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summaryIntPx', VarColIndex=2, YLimMaxSet = 110,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCTCCM_Glass_Class2Int'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryIntPx(:,3), VarNaming = 'Small Circular', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='yellows', VarColName='FA_summaryIntPx', VarColIndex=3, YLimMaxSet = 110,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCTCCM_Glass_Class3Int'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryIntPx(:,4), VarNaming = 'Large Circular', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='purples', VarColName='FA_summaryIntPx', VarColIndex=4, YLimMaxSet = 110,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCTCCM_Glass_Class4Int'), 'svg', Width=16, Height=12);
%% 
plotMetricScatter(T_Vin, T_Vin.FA_summaryFRETPx(:,1), VarNaming = 'Linear', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=15, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summaryFRETPx', VarColIndex=1, YLimMaxSet = 35,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCTCCM_Glass_Class1FRET'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFRETPx(:,2), VarNaming = 'Elliptical', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=15, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summaryFRETPx', VarColIndex=2, YLimMaxSet = 35,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCTCCM_Glass_Class2FRET'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFRETPx(:,3), VarNaming = 'Small Circular', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=15, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='yellows', VarColName='FA_summaryFRETPx', VarColIndex=3, YLimMaxSet = 35,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCTCCM_Glass_Class3FRET'), 'svg', Width=16, Height=12);
plotMetricScatter(T_Vin, T_Vin.FA_summaryFRETPx(:,4), VarNaming = 'Large Circular', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=15, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='purples', VarColName='FA_summaryFRETPx', VarColIndex=4, YLimMaxSet = 35,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCTCCM_Glass_Class4FRET'), 'svg', Width=16, Height=12);
%% Figure S22: Intensity, FRET between classification
plotMetricScatter(longTblInt, 'FA_summaryIntPx', excludePatterns={'CCM'}, includePatterns={'1kPa'}, YLimMinSet=0, YLimMaxSet = 105,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='1kPa',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCT_1kPa_Intcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblInt, 'FA_summaryIntPx', excludePatterns={'CT'}, includePatterns={'1kPa'}, YLimMinSet=0, YLimMaxSet = 105,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='1kPa',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCCM_1kPa_Intcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblInt, 'FA_summaryIntPx', excludePatterns={'CCM'}, includePatterns={'10kPa'}, YLimMinSet=0, YLimMaxSet = 105,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='10kPa',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCT_10kPa_Intcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblInt, 'FA_summaryIntPx', excludePatterns={'CT'}, includePatterns={'10kPa'}, YLimMinSet=0, YLimMaxSet = 105,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='10kPa',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCCM_10kPa_Intcat'), 'svg', Width=16, Height=12);
%%
plotMetricScatter(longTblFRET, 'FA_summaryFRETPx', excludePatterns={'CCM'}, includePatterns={'1kPa'}, YLimMinSet=15, YLimMaxSet = 37,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='1kPa',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCT_1kPa_FRETcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblFRET, 'FA_summaryFRETPx', excludePatterns={'CT'}, includePatterns={'1kPa'}, YLimMinSet=15, YLimMaxSet = 37,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='1kPa',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCCM_1kPa_FRETcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblFRET, 'FA_summaryFRETPx', excludePatterns={'CCM'}, includePatterns={'10kPa'}, YLimMinSet=15, YLimMaxSet = 37,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='10kPa',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCT_10kPa_FRETcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblFRET, 'FA_summaryFRETPx', excludePatterns={'CT'}, includePatterns={'10kPa'}, YLimMinSet=15, YLimMaxSet = 37,...
    XCat='fractionCat', levels=Input.CatNames, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='10kPa',...
    CondColors=Input.condColorsClass, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig6_SI_VinCCM_10kPa_FRETcat'), 'svg', Width=16, Height=12);
%% SI - Extra:
%% Differences plots
% plotMetricScatter(T_Vin, 'Av_TensionPerDens', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0,  ...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors,...
%     ConditionOrder = Input.CondOrderGel);
% plotMetricScatter(T_diffVinkPa, 'FA_Count', levels={'Diff_10kPav1kPa'}, excludePatterns={}, YLimMinSet=0,  ...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors, statTest2 ='wlm');
% plotMetricScatter(T_diffVinkPa_Rev, 'FRETav', levels={'Diff_1kPav10kPa'}, excludePatterns={}, YLimMinSet=0,  ...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors, statTest2 ='wlm');
% plotMetricScatter(T_diffVinkPa, 'Av_AverageIntensity', levels={'Diff_10kPav1kPa'}, excludePatterns={}, YLimMinSet=0,  ...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors, statTest2 ='wlm');
% plotMetricScatter(T_diffVinkPa2, 'FA_Count', levels={'Diff_CCMvCT'}, excludePatterns={}, YLimMinSet=0,  ...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors, statTest2 ='wlm',...
%     XCat='diseaseCat', compareCat='extrainfoCat');
% plotMetricScatter(T_diffVinkPa2, 'FRETav', levels={'Diff_CCMvCT'}, excludePatterns={}, ...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors, statTest2 ='wlm',...
%     XCat='diseaseCat', compareCat='extrainfoCat');
% plotMetricScatter(T_diffVinkPa2, 'Av_AverageIntensity', levels={'Diff_CCMvCT'}, excludePatterns={},...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors, statTest2 ='wlm',...
%     XCat='diseaseCat', compareCat='extrainfoCat');


% ====================================================================
%% ----------------- FIG_7: VECad on Gel  -----------------
% ====================================================================
%% 12 Pannel (4 sets)
PaperImg_ExportIntFRETimage(T_VECad, Input, '20260429_VECadTS_FL_CT_1kPa_si72', 15, ...
    Class=true, Skel=true, ThickMax=11, IntLow = 0, IntHigh = 99.99, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_VECad, Input, '20260427_VECadTS_FL_CT_10kPa_si72', 12, ...
    Class=true, Skel=true, ThickMax=11, IntLow = 0, IntHigh = 99.5, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_VECad, Input, '20260429_VECadTS_FL_CCM_1kPa_si72', 9, ...
    Class=true, Skel=true, ThickMax=11, IntLow = 0, IntHigh = 99.80, ExportSVG = true)
PaperImg_ExportIntFRETimage(T_VECad, Input, '20260427_VECadTS_FL_CCM_10kPa_si72', 6, ...
    Class=false, Skel=false, ThickMax=11, IntLow = 0, IntHigh = 99.9, ExportSVG = true)
%% CellShape, Thickness(skel), Intensity, FRET
plotMetricScatter(T_VECad, 'CellAv_MinMajRatio', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, Unit='', ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors,...
    ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_VECadCTCCM_Gel_CellShape'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, 'Thickness_Skeleton_Av', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, Unit='(µm)', ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors,...
    ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_VECadCTCCM_Gel_Thick'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, 'Av_AverageIntensity', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    Unit='', ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove,...
    CondColors=Input.condColors, ThirdCatFilter={'TS'}, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_VECadCTCCM_Gel_AvInt'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, 'FRETav', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=15, Unit='(%) - Apparent', ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors,...
    ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_VECadCTCCM_Gel_FRETav'), 'svg', Width=16, Height=12);
%% Histogram
Plot_ShapeClassStackedHistogram(T_VECad, Levels={'1kPa', '10kPa'}, excludePatterns={'TL'}, legendNames = Input.CatNamesAJ)
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_VECadCTCCM_Gel_ClassHist'), 'svg');

%% SI
%% Figure S24: Cell Size
plotMetricScatter(T_VECad, 'CellAvUmSize', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, Unit='(µm^2)', ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, ContrastRemove = Input.ContrastRemove, CondColors=Input.condColors,...
    ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCTCCM_Gel_CellSize'), 'svg', Width=16, Height=12);
%% Figure S26: Fraction, intensity, FRET per classification
plotMetricScatter(T_VECad, T_VECad.FA_summaryFrac(:,1), VarNaming = 'Reticular', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summaryFrac', VarColIndex=1, YLimMaxSet = 1.2,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCTCCM_Gel_Class1Frac'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, T_VECad.FA_summaryFrac(:,2), VarNaming = 'Linear', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summaryFrac', VarColIndex=2, YLimMaxSet = 1.2,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCTCCM_Gel_Class2Frac'), 'svg', Width=16, Height=12);
%%
plotMetricScatter(T_VECad, T_VECad.FA_summarySkelInt(:,1), VarNaming = 'Reticular', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summarySkelInt', VarColIndex=1, YLimMaxSet = 200,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCTCCM_Gel_Class1Int'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, T_VECad.FA_summarySkelInt(:,2), VarNaming = 'Linear', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=0, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summarySkelInt', VarColIndex=2, YLimMaxSet = 200,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCTCCM_Gel_Class2Int'), 'svg', Width=16, Height=12);
%% 
plotMetricScatter(T_VECad, T_VECad.FA_summarySkelFRET(:,1), VarNaming = 'Reticular', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=15, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='bluepurples', VarColName='FA_summarySkelFRET', VarColIndex=1, YLimMaxSet = 38,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCTCCM_Gel_Class1FRET'), 'svg', Width=16, Height=12);
plotMetricScatter(T_VECad, T_VECad.FA_summarySkelFRET(:,2), VarNaming = 'Linear', levels={'1kPa', '10kPa'}, excludePatterns={}, YLimMinSet=15, ...
    ContrastExcludePatterns = Input.ContrastExcludePatterns, paletteName='greens', VarColName='FA_summarySkelFRET', VarColIndex=2, YLimMaxSet = 38,...
    ContrastRemove = Input.ContrastRemove, ConditionOrder = Input.CondOrderGel);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCTCCM_Gel_Class2FRET'), 'svg', Width=16, Height=12);
%% Figure S27: Intensity, FRET between classification
plotMetricScatter(longTblIntVE, 'FA_summarySkelInt', excludePatterns={'CCM'}, includePatterns={'1kPa'}, YLimMinSet=0, YLimMaxSet = 160,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='1kPa',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCT_1kPa_Intcat'), 'svg', Width=16, Height=12);;
plotMetricScatter(longTblIntVE, 'FA_summarySkelInt', excludePatterns={'CT'}, includePatterns={'1kPa'}, YLimMinSet=0, YLimMaxSet = 160,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='1kPa',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCCM_1kPa_Intcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblIntVE, 'FA_summarySkelInt', excludePatterns={'CCM'}, includePatterns={'10kPa'}, YLimMinSet=0, YLimMaxSet = 160,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='10kPa',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCT_10kPa_Intcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblIntVE, 'FA_summarySkelInt', excludePatterns={'CT'}, includePatterns={'10kPa'}, YLimMinSet=0, YLimMaxSet = 160,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='10kPa',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCCM_10kPa_Intcat'), 'svg', Width=16, Height=12);
%%
plotMetricScatter(longTblFRETVE, 'FA_summarySkelFRET', excludePatterns={'CCM'}, includePatterns={'1kPa'}, YLimMinSet=15, YLimMaxSet = 38,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='1kPa',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCT_1kPa_FRETcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblFRETVE, 'FA_summarySkelFRET', excludePatterns={'CT'}, includePatterns={'1kPa'}, YLimMinSet=15, YLimMaxSet = 38,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='1kPa',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCCM_1kPa_FRETcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblFRETVE, 'FA_summarySkelFRET', excludePatterns={'CCM'}, includePatterns={'10kPa'}, YLimMinSet=15, YLimMaxSet = 38,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='10kPa',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCT_10kPa_FRETcat'), 'svg', Width=16, Height=12);
plotMetricScatter(longTblFRETVE, 'FA_summarySkelFRET', excludePatterns={'CT'}, includePatterns={'10kPa'}, YLimMinSet=15, YLimMaxSet = 38,...
    XCat='fractionCat', levels=Input.CatNamesAJ, compareCat='diseaseCat', SplitCat='extrainfoCat', SplitCatLevel='10kPa',...
    CondColors=Input.condColorsClassAJ, statTest2='rm');
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig7_SI_VECadCCM_10kPa_FRETcat'), 'svg', Width=16, Height=12);
%% SI - Extra
%% Differences plots
% plotMetricScatter(T_diffVEkPa_Rev, 'Thickness_Skeleton_Av', levels={'Diff_1kPav10kPa'}, excludePatterns={}, YLimMinSet=0,  ...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors, statTest2 ='wlm');
% plotMetricScatter(T_diffVEkPa_Rev, 'FRETav', levels={'Diff_1kPav10kPa'}, excludePatterns={}, YLimMinSet=0,  ...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors, statTest2 ='wlm');
% plotMetricScatter(T_diffVEkPa, 'Av_AverageIntensity', levels={'Diff_10kPav1kPa'}, excludePatterns={}, YLimMinSet=0,  ...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors, statTest2 ='wlm');
% plotMetricScatter(T_diffVEkPa2, 'FA_Count', levels={'Diff_CCMvCT'}, excludePatterns={}, YLimMinSet=0,  ...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors, statTest2 ='wlm',...
%     XCat='diseaseCat', compareCat='extrainfoCat');
% plotMetricScatter(T_diffVEkPa2, 'FRETav', levels={'Diff_CCMvCT'}, excludePatterns={}, ...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors, statTest2 ='wlm',...
%     XCat='diseaseCat', compareCat='extrainfoCat');
% plotMetricScatter(T_diffVEkPa2, 'Av_AverageIntensity', levels={'Diff_CCMvCT'}, excludePatterns={},...
%     ContrastExcludePatterns = Input.ContrastExcludePatterns, CondColors=Input.condColors, statTest2 ='wlm',...
%     XCat='diseaseCat', compareCat='extrainfoCat');


% ====================================================================
%% ----------------- FIG_8: Vin&VECad heatmaps  -----------------
% ====================================================================
%% First combine the 2 results tables into one table to generate one heatmap
% Rename columns 4+ in both tables by adding 'FA_' prefix
comparisonResultsFA.Properties.VariableNames(4:end) = strcat('FA_', comparisonResultsFA.Properties.VariableNames(4:end));
% Extract columns 4+ from table2 and append to table1
allcomparisonResults = [comparisonResultsFA, comparisonResultsAJ(:, 4:end)];
%% Plot into one heatmap
Plot_ConditionsCohensD_heatmap(allcomparisonResults, Input.CompareVariables);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig8_Heat'), 'svg', Width=20, Height=14);
%saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig8_Heat'), 'svg', Width=25, Height=16);
%%
Plot_ConditionsCohensD_heatmap(comparisonResultsFA, Input.CompareVariablesFA, Levels = {'CCM_rest', 'CCM_R1', 'CCM_R2'}, CLim=[-3 3]);
saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig8_Heat_FAROCK'), 'svg', Width=15, Height=9);
Plot_ConditionsCohensD_heatmap(comparisonResultsAJ, Input.CompareVariablesAJ, Levels = {'CCM_rest', 'CCM_R1', 'CCM_R2'}, CLim=[-3 3]);
saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig8_Heat_AJROCK'), 'svg', Width=15, Height=9);
%%
Plot_ConditionsCohensD_heatmap(comparisonResultsFA, Input.CompareVariablesFA, Levels = {'CT_1kPa', 'CT_10kPa', 'CCM_1kPa', 'CCM_10kPa'});
saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig8_Heat_FAGel'), 'svg', Width=15, Height=9);
Plot_ConditionsCohensD_heatmap(comparisonResultsAJ, Input.CompareVariablesAJ, Levels = {'CT_1kPa', 'CT_10kPa', 'CCM_1kPa', 'CCM_10kPa'});
saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig8_Heat_AJGel'), 'svg', Width=15, Height=9);
%%
Plot_ConditionsCohensD_heatmap(comparisonResultsFA, Input.CompareVariablesFA, Levels = {'CT_1kPa', 'CCM_1kPa', 'CT_10kPa', 'CCM_10kPa'});
saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig8_Heat_FAGelOrder'), 'svg', Width=15, Height=9);
Plot_ConditionsCohensD_heatmap(comparisonResultsAJ, Input.CompareVariablesAJ, Levels = {'CT_1kPa', 'CCM_1kPa', 'CT_10kPa', 'CCM_10kPa'});
saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig8_Heat_AJGelOrder'), 'svg', Width=15, Height=9);

%% Vin
% Plot_ConditionsCohensD_heatmap(comparisonResultsFA, Input.CompareVariablesFA, TitleInfo='— Vin');
% %saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig8_VinHeat'), 'svg', Width=20, Height=14);
% % Plot_ConditionsCohensD_heatmap(allcomparisonResultsFA.CatAll, Input.CompareVariablesFA, TitleInfo='— All');
% % Plot_ConditionsCohensD_heatmap(allcomparisonResultsFA.TCat1,  Input.CompareVariablesFA, TitleInfo='— Cat1');
% % Plot_ConditionsCohensD_heatmap(allcomparisonResultsFA.TCat2,  Input.CompareVariablesFA, TitleInfo='— Cat2');
% % Plot_ConditionsCohensD_heatmap(allcomparisonResultsFA.TCat3,  Input.CompareVariablesFA, TitleInfo='— Cat3');
% % Plot_ConditionsCohensD_heatmap(allcomparisonResultsFA.TCat4,  Input.CompareVariablesFA, TitleInfo='— Cat4');
%% VECad
% Plot_ConditionsCohensD_heatmap(comparisonResultsAJ, Input.CompareVariablesAJ, TitleInfo='— VECad');
% saveVectorFigure(gcf, fullfile(Input.OutPath, 'Fig8_VECadHeat'), 'svg', Width=20, Height=14);
% % Plot_ConditionsCohensD_heatmap(allcomparisonResultsAJ.CatAll, Input.CompareVariablesAJ, TitleInfo='— All');
% % Plot_ConditionsCohensD_heatmap(allcomparisonResultsAJ.TCat1,  Input.CompareVariablesAJ, TitleInfo='— Cat1');
% % Plot_ConditionsCohensD_heatmap(allcomparisonResultsAJ.TCat2,  Input.CompareVariablesAJ, TitleInfo='— Cat2');

%% IntFRET
% % Assume you already have 512x512 arrays:
% % mask (binary), ch1 (uint8 or double), ch2 (uint8 or double)
% 
% mask = Cell1_DD_ch00_FASegment;
% ch1 = Cell1_DD_ch00;
% ch2 = Cell1_DA_ch01;
% FRET_ColorbarLim = 1;
% 
% T = IntFRET(mask, ch1, ch2);
% 
% % Access the FRET image:
% F = T.FRETimage{1};
% %figure; imagesc(F); axis image off; colorbar; title('FRET ratio (ch2/ch1)'); xlim([0 user.FRET_ColorbarLim]);
% figure;
% imagesc(F, [0 1]);   % Force the displayed range to 0–1
% axis image off;
% colormap(jet);       % Or 'turbo', 'parula', etc.
% colorbar;
% caxis([0 1]);        % Ensures the colorbar also spans 0–1
% title('FRET ratio (0–1 scaled display)');
% 
% % Average FRET per object (Nx1):
% F_obj = T.FRETperObject{1};
% 
% % Whole-image average FRET:
% F_mean = T.FRETav;

%% PCA
% Input.VinPCAColumns = {'Area','Ecc','MajAx','MinAx','MinMajRatio','NND3','FADens_r5'};
% Input.VinPCAColumns = {'Area','Ecc','MajAx','MinMajRatio','NND3','IntensityPerPixel', 'MaxInt'};
% % [X, varNames, coeff, score, explained, clusterIdx, silh] = ...
% %     runPCA_EccResults_3D(T_Vin, Input);
% % RESULT: 3D PCA necessary (PC3 = ~13%)
% %%
% [X, varNames, coeffTbl, scoreTbl, explainedTbl, usedMask, groupLabels] = ...
%     runPCA_FOV_Aggregated(T_Vin, Input, 'median', true);
% %%
% [X, varNames, coeffTbl, scoreTbl, explainedTbl, usedMask, groupLabels] = ...
% runPCA_FOV_Aggregated_Filtered(T_Vin, Input, "mean", false, ["CT-rest", "CCM-rest"]);



