function Save2Excel(outputMetrics, pointers, countLoop, flags)
%SAVE2EXCEL Summary of this function goes here
%   Detailed explanation goes here

%% Create excel file

filename = [outputMetrics.outputPath filesep 'FDA_Vinculin_' pointers.DirName '.xlsx'];
    
if countLoop == 1 % only do this once

    SampleName = {'Name:' ; 'Path:'};
    SampleFolder = {string(pointers.DirName) ; string(pointers.directory)};

    T = table({'test--------' ; 'one'}, {'test--------' ; 'two'});
    writetable(T, filename,'Sheet', 'Pooled cells','Range','A1', 'WriteVariableNames', false)

    T = table(SampleName, SampleFolder);
    writetable(T, filename,'Sheet', 'Pooled cells','Range','A1', 'WriteVariableNames', false, "AutoFitWidth", false);%, "WriteMode","append")

else
end

%% write average FRET

% create table with Headers
T = table({'Average FRET FA'}, {'Average FRET Cytosol'}, {'Average FRET BG'}, {'Cell area (µm²)'},...
    {'Cell intensity'}, {'Sum FA intensity'}, {'Cytosol intensity'}, {'Average FA Int per Cell'});

% write table to excel
cellCoords = ['B', num2str(4)];
writetable(T,filename,'Sheet', 'Pooled cells','Range', cellCoords, 'WriteVariableNames', false, "AutoFitWidth", false)

% parse
FileName = pointers.FileName;
FRETav = outputMetrics.FRETav;
FRETav_cyto = outputMetrics.FRETav_cyto;
FRETav_bg = outputMetrics.FRETav_bg;
cellAreaUm2 = outputMetrics.cellAreaUm2;
cellIntensity = outputMetrics.cellIntensity;
summedFAIntensity = outputMetrics.totalFAintensity;
cytosolIntensity = outputMetrics.cytosolIntensity;
meanOfAllFAintensity = outputMetrics.meanOfTotalFAintensity;

% create table with values
TC = {convertCharsToStrings(FileName) str2num(FRETav) str2num(FRETav_cyto) str2num(FRETav_bg) cellAreaUm2 cellIntensity summedFAIntensity cytosolIntensity, meanOfAllFAintensity};

T = cell2table(TC);

% write table to excel
cellCoords = ['A', num2str(countLoop + 4)];
writetable(T,filename,'Sheet', 'Pooled cells','Range', cellCoords, 'WriteVariableNames', false, "AutoFitWidth", false)

%% write metrics
AverageFRETperFA = outputMetrics.AverageFRETperFA;
AreaInMicroMeterSq = outputMetrics.AreaInMicroMeterSq;
AreaInPixels = outputMetrics.AreaInPixelsFA;
TotalIntensity = outputMetrics.TotalIntensity;
AverageIntensityPerPixel = outputMetrics.AverageIntensity;
FA_index = (1:1:length(AverageFRETperFA))';
AverageIntensityPerMicroMeterSq = outputMetrics.IntensityPerMicroMeterSq;

T = table(FA_index, AverageFRETperFA, AreaInMicroMeterSq, AreaInPixels, TotalIntensity, AverageIntensityPerPixel, AverageIntensityPerMicroMeterSq);
sheetName = pointers.FileName;
writetable(T,filename,'Sheet', sheetName,'Range','A1')

%% --- Write histogram data for this cell ---
if flags.flagPhasorHistogram
histTabName = [pointers.FileName ' Hist'];  % e.g., "Cell1 Hist"
Histogram_BinCenters = outputMetrics.Histogram_BinCenters;
Histogram_Counts = outputMetrics.Histogram_Counts;
Histogram_NormalizedCounts = outputMetrics.Histogram_NormalizedCounts;
% Create table
T_Hist = table(Histogram_BinCenters', ...
               Histogram_Counts', ...
               Histogram_NormalizedCounts', ...
               'VariableNames', {'BinCenter', 'Counts', 'NormalizedCounts'});

% Write table to new sheet for this cell
writetable(T_Hist, filename, 'Sheet', histTabName, 'Range', 'A1');
end
end

