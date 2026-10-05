function LineAnalysis2Excel(data, FA, pointers)
%LINEANALYSIS3EXCEL Save FRET and Intensity line analysis for one cell
% Each call saves one cell, creating sheets named after the cell

outPath = data.outputPath;

% Always define timestamp string (used for Date column in Excel)
t = datetime;
t.Format = 'yyMMdd-HHmm';
t = convertStringsToChars(string(t));

% Look for an existing Line_Analysis_ file in the output folder
existingFile = dir(fullfile(outPath, 'Line_Analysis_*.xlsx'));

if ~isempty(existingFile)
    % Reuse the first match (or pick newest one if desired)
    filename = fullfile(outPath, existingFile(1).name);
else
    % If none exists, create a new one with current timestamp
    filenamePath = [outPath filesep 'Line_Analysis_' t];
    filename = [filenamePath '.xlsx'];
end


startCell = 'A2';

% Dynamic sheet names based on cell
sheetFRET = [pointers.FileName '_FRET'];
sheetInt = [pointers.FileName '_Intensity'];

%% unique ID array
Perimeter_ID = FA.Perimeter;
FA_ID = (1:1:length(Perimeter_ID))'; 
Cell_ID = [data.SampleDate '_' data.name];
Cell_ID = repmat({Cell_ID}, length(FA_ID), 1);
Date = repmat({t}, length(FA_ID), 1);
Raw = repmat({'Raw'}, length(FA_ID), 1); 
Interpolated = repmat({'Interpolated'}, length(FA_ID), 1); 
Normalized = repmat({'Normalized'}, length(FA_ID), 1); 
%% FA-Line-FRET
Len = cellfun(@length, FA.intra.lineFaIdxList);
[lineLengthFApix] = FillLines(Len', 0, FA.intra.interp.FRET);

var.FA1 = Cell_ID;
var.FA2 = FA_ID;
var.FA3 = Perimeter_ID;
var.FA4 = lineLengthFApix;
var.Raw = Raw;
var.Interpolated = Interpolated;
var.Normalized = Normalized; 
var.avg = FA.intra.FRETpadded;
%var.header1 = {'FRETLineUntouched'};
var.header1 = 1:50; 
var.header2 = {'Interpolated'};
var.interp = FA.intra.interp.FRET;
var.interpAverage = FA.intra.interp.FRETaverage;
var.header3 = {'AverageInterpolatedLine'};
var.norm = FA.intra.interp.normFRET;
var.normAverage = FA.intra.interp.normFRETaverage;
var.header4 = {'Normalized'};
var.header5 = {'AverageNormalizedLine'};

Write2excelLines(var, filename, sheetFRET, startCell);
var = [];

%% FA-Line-Intensity
var.FA1 = Cell_ID;
var.FA2 = FA_ID;
var.FA3 = Perimeter_ID;
var.FA4 = lineLengthFApix;
var.Raw = Raw;
var.Interpolated = Interpolated;
var.Normalized = Normalized; 
var.avg = FA.intra.Intensitypadded;
%var.header1 = {'IntensityLineUntouched'};
var.header1 = 1:50; 
var.header2 = {'Interpolated'};
var.interp = FA.intra.interp.intensity;
var.interpAverage = FA.intra.interp.intensityaverage;
var.header3 = {'AverageInterpolatedLine'};
var.norm = FA.intra.interp.normIntensity;
var.normAverage = FA.intra.interp.normIntensityaverage;
var.header4 = {'Normalized'};
var.header5 = {'AverageNormalizedLine'};

Write2excelLines(var, filename, sheetInt, startCell);
var = [];

end


%% --- Subfunctions ---
function [output1] = FillLines(target, filler, targetLength)
m = 1;
for k = 1:length(target)  % <- loop over number of FAs, NOT size(targetLength,1)
    if sum(targetLength(k, :)) == 0
        output1(k,1) = filler;
    else
        output1(k,1) = target(m, 1);
        m = m + 1;
    end
end
end


function Write2excelLines(var, filename, sheetName, startRange)

T1raw = table(var.FA1, var.FA2, var.FA3, var.FA4, var.Raw, ...
    'VariableNames', {'SampleName','FA_ID','PerimeterTag','LineLength', 'LineType'});
T1interpolated = table(var.FA1, var.FA2, var.FA3, var.FA4, var.Interpolated, ...
    'VariableNames', {'SampleName','FA_ID','PerimeterTag','LineLength', 'LineType'});
T1normalized = table(var.FA1, var.FA2, var.FA3, var.FA4, var.Normalized, ...
    'VariableNames', {'SampleName','FA_ID','PerimeterTag','LineLength', 'LineType'});
T2 = table(var.avg);
T3 = table(var.header1);
T4 = table(var.header2);
T5 = table(var.interp);
T6 = table(var.interpAverage);
T7 = table(var.header3);
T8 = table(var.norm);
T9 = table(var.normAverage);
T10 = table(var.header4);
T11 = table(var.header5);

% Start writing
% First we write the first 4 columns
writetable(T1raw, filename, 'Sheet', sheetName, 'Range', startRange)

alphabet = string(('A':'Z').').';
% Then we write the untouched FRET padded data
startCell_tmp2 = startRange;
startCell_tmp2(1) =  alphabet(size(T1raw, 2)+1);
startCell_tmp2(2) = num2str(str2double(startCell_tmp2(2)) + 1); 
writetable(T2, filename, 'Sheet', sheetName, 'Range', startCell_tmp2, 'WriteVariableNames', false)
% Then we put the title of untouched data above it
startCell_tmp3 = startRange;
startCell_tmp3(1) = alphabet(size(T1raw, 2)+1);
writetable(T3, filename, 'Sheet', sheetName, 'Range', startCell_tmp3, 'WriteVariableNames', false)
%Then we write the first 4 column again
startCell_tmp1b = startRange;
startCell_tmp1b = [startCell_tmp1b(1) num2str(str2double(startCell_tmp1b(2:end)) + size(var.FA1, 1) + 4)];
writetable(T1interpolated, filename, 'Sheet', sheetName, 'Range', startCell_tmp1b, 'WriteVariableNames', false)
%Then we first put the title of the interpolated data this time 
startCell_tmp4 = startCell_tmp1b;
startCell_tmp4 = [startCell_tmp4(1) num2str(str2double(startCell_tmp4(2:end)) - 1)];
startCell_tmp4(1) = alphabet(size(T1interpolated, 2)+1);
writetable(T4, filename, 'Sheet', sheetName, 'Range', startCell_tmp4, 'WriteVariableNames', false)
% Then we finally put the interpolated data itsel 
startCell_tmp5 = startCell_tmp4;
startCell_tmp5 = [startCell_tmp5(1) num2str(str2double(startCell_tmp5(2:end)) + 1)];
writetable(T5, filename, 'Sheet', sheetName, 'Range', startCell_tmp5, 'WriteVariableNames', false)
% Then we put the average interpolated data
startCell_tmp6 = startCell_tmp5;
startCell_tmp6 = [startCell_tmp6(1) num2str(str2double(startCell_tmp6(2:end)) + size(var.FA1, 1) + 1)];
writetable(T6, filename, 'Sheet', sheetName, 'Range', startCell_tmp6, 'WriteVariableNames', false)
% Then we put the title of the average interpolated data above 
startCell_tmp7 = startCell_tmp6;
startCell_tmp7 = [startCell_tmp7(1) num2str(str2double(startCell_tmp7(2:end)) - 1)];
writetable(T7, filename, 'Sheet', sheetName, 'Range', startCell_tmp7, 'WriteVariableNames', false)
% Then we put the first four column one last time for the normalized data
startCell_tmp1c = startCell_tmp1b;
startCell_tmp1c = ['A' num2str(str2double(startCell_tmp6(2:end)) + 4)];
writetable(T1normalized, filename, 'Sheet', sheetName, 'Range', startCell_tmp1c, 'WriteVariableNames', false)
% 
startCell_tmp10 = startCell_tmp1c;
startCell_tmp10 = [startCell_tmp10(1) num2str(str2double(startCell_tmp10(2:end)) - 1)];
startCell_tmp10(1) = alphabet(size(T1normalized, 2)+1);
writetable(T10, filename, 'Sheet', sheetName, 'Range', startCell_tmp10, 'WriteVariableNames', false)

startCell_tmp8 = startCell_tmp10;
startCell_tmp8 = [startCell_tmp8(1) num2str(str2double(startCell_tmp8(2:end)) + 1)];
writetable(T8, filename, 'Sheet', sheetName, 'Range', startCell_tmp8, 'WriteVariableNames', false)

startCell_tmp9 = startCell_tmp8;
startCell_tmp9 = [startCell_tmp9(1) num2str(str2double(startCell_tmp9(2:end)) + size(var.FA1, 1) + 1)];
writetable(T9, filename, 'Sheet', sheetName, 'Range', startCell_tmp9, 'WriteVariableNames', false)

startCell_tmp11 = [startCell_tmp9(1) num2str(str2double(startCell_tmp9(2:end)) - 1)];
writetable(T11, filename, 'Sheet', sheetName, 'Range', startCell_tmp11, 'WriteVariableNames', false)

end
