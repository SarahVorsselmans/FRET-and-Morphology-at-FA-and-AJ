function T = Order_TableByRules(T, Info)
% Sort_Table_ByRules
% -------------------------------------------------------------
% Applies categorical ordering rules from Info to table T and
% sorts the table by:
%   spring -> protein -> sampletype -> disease -> extrainfo -> date
%
% Expected naming conventions:
%   SampleFolder example:
%   20260125_VinTS_FL_CT_R2_si72
%   20260125_VECadTL_FL_CCM_si72
%
% Protein extraction rules:
%   - "Vin" → protein = Vin
%   - "VECad" → protein = VECad
%
% REQUIRED Info fields:
%   springOrder
%   proteinOrder    (e.g. {'Vin','VECad'})
%   sampletypeOrder (e.g. {'TS','TL'})
%   diseaseOrder
%   extrainfoOrder
% -------------------------------------------------------------

%% ---------- EXPERIMENT ----------
T.expCat = categorical(T.DataFolder);

%% ---------- DATE ----------
T.dateCat = categorical(T.date);

%% ---------- SPRING ----------
T.springCat = categorical(T.spring, Info.springOrder, 'Ordinal', true);

%% ---------- PROTEIN (always Vin or VECad) ----------
proteinStr = strings(height(T),1);

isVin = contains(T.SampleFolder, "_Vin", "IgnoreCase", true);
isVECad  = contains(T.SampleFolder, "_VECad",  "IgnoreCase", true);

proteinStr(isVin) = "Vin";
proteinStr(isVECad)  = "VECad";

% Convert to categorical with your defined order
T.proteinCat = categorical(proteinStr, Info.proteinOrder, 'Ordinal', true);

%% ---------- SAMPLETYPE ----------
% sampletype is something like "VinTS", "VECadTL", etc.
% But we only need TS vs TL

isTS = contains(T.sampletype, "TS", "IgnoreCase", true);

stype = repmat("TL", height(T),1);
stype(isTS) = "TS";

T.sampletypeCat = categorical(stype, Info.sampletypeOrder, 'Ordinal', true);

%% ---------- DISEASE ----------
T.diseaseCat = categorical(T.disease, Info.diseaseOrder, 'Ordinal', true);

%% ---------- EXTRAINFOS ----------
% Collapse multiple variants into clean categories
exinfo = repmat("rest", height(T),1);

exinfo(contains(T.extrainfo,"R1"))    = "R1";
exinfo(contains(T.extrainfo,"R2"))    = "R2";
exinfo(contains(T.extrainfo,"10kPa")) = "10kPa";
exinfo(contains(T.extrainfo,"1kPa"))  = "1kPa";

T.extrainfoCat = categorical(exinfo, Info.extrainfoOrder, 'Ordinal', true);

% %% ---------- LIVE vs FIXED ----------
% Collapse multiple variants into clean categories
% liveinfo = repmat("Fixed", height(T),1);
% 
% liveinfo(contains(T.extrainfo,"live"))  = "Live";
% 
% T.LiveFixedCat = categorical(liveinfo, Info.LiveFixedOrder, 'Ordinal', true);

%% ---------- xCat (SampleFolder stable order) ----------
T.xCat = categorical(T.SampleFolder, unique(T.SampleFolder,'stable'), 'Ordinal', true);

%% ---------- Build sorting table ----------
sortTable = table( ...
    T.springCat, ...
    T.proteinCat, ...
    T.sampletypeCat, ...
    T.diseaseCat, ...
    T.extrainfoCat, ...
    categorical(T.date) ...
);

[~, sortIdx] = sortrows(sortTable);
T = T(sortIdx, :);

%% ---------- Fix xCat order according to sorted T ----------
T.xCat = reordercats(T.xCat, cellstr(unique(T.xCat, 'stable')));

end