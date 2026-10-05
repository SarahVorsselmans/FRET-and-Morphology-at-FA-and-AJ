function [T, BaselineSummary] = Norm_ToCTBaseline(T, metricField, targetCT, combinePureCTAcrossBaseExtra)
%NORMALIZETOCTBASELINE Normalize a metric to CT baseline per DataFolder (optionally combining across baseExtra).
%   [T, BaselineSummary] = NormalizeToCTBaseline(T)
%   [T, BaselineSummary] = NormalizeToCTBaseline(T, metricField, targetCT, combinePureCTAcrossBaseExtra)
%
% Inputs:
%   T           - table with variables:
%                 DataFolder, extrainfo, diseaseCat, sampletypeCat, springCat,
%                 and the metric to normalize (default: 'FRETav').
%   metricField - (char/string) metric field name in T to normalize (default: 'FRETav').
%   targetCT    - scalar to scale CT mean to (default: 1).
%   combinePureCTAcrossBaseExtra - logical (default: true).
%                 If true, baseline is per (DataFolder, sampletypeCat, springCat), combining all pure CT
%                 regardless of baseExtra (e.g., si72 and glass are pooled for the baseline).
%                 If false, baseline also includes baseExtra (original behavior).
%
% Outputs:
%   T               - table with added columns:
%                     T.(metricField + "ctMeanNorm")  (per-group CT baseline)
%                     T.(metricField + "Norm")        (normalized metric scaled to targetCT)
%   BaselineSummary - summary of baselines used per group. If combining across baseExtra,
%                     the summary includes the count of pure CT rows pooled from all baseExtra.
%
% Notes:
%   - "Pure CT" means diseaseCat indicates CT and extrainfo does not contain R# or kPa tokens.
%   - Grouping always includes (DataFolder, sampletypeCat, springCat).
%     baseExtra is used for diagnostics or (optionally) as part of the grouping for baseline.

    % ---- Defaults ----
    if nargin < 2 || isempty(metricField), metricField = 'FRETav'; end
    if nargin < 3 || isempty(targetCT),    targetCT    = 1;        end
    if nargin < 4 || isempty(combinePureCTAcrossBaseExtra), combinePureCTAcrossBaseExtra = true; end

    metricField = string(metricField);

    % ---- Dynamic variable names ----
    ctMeanVar = matlab.lang.makeValidName(metricField + "ctMeanNorm");
    normVar   = matlab.lang.makeValidName(metricField + "Norm");

    % ---- Basic validation ----
    mustHave = {'DataFolder','extrainfo','diseaseCat','sampletypeCat','springCat', metricField};
    for i = 1:numel(mustHave)
        if ~ismember(mustHave{i}, T.Properties.VariableNames)
            error('NormalizeToCTBaseline:MissingVar', ...
                  'Variable "%s" not found in T.', mustHave{i});
        end
    end
    if ~isnumeric(T.(metricField))
        error('NormalizeToCTBaseline:ValueNotNumeric', ...
              'The metric "%s" must be numeric.', metricField);
    end
    if ~iscategorical(T.diseaseCat)
        warning('NormalizeToCTBaseline:Type', 'T.diseaseCat is expected categorical; converting.');
        T.diseaseCat = categorical(T.diseaseCat);
    end
    if ~iscategorical(T.sampletypeCat)
        warning('NormalizeToCTBaseline:Type', 'T.sampletypeCat is expected categorical; converting.');
        T.sampletypeCat = categorical(T.sampletypeCat);
    end
    if ~iscategorical(T.springCat)
        warning('NormalizeToCTBaseline:Type', 'T.springCat is expected categorical; converting.');
        T.springCat = categorical(T.springCat);
    end

    % ---- Keys ----
    dataFolderStr = string(T.DataFolder);
    rawExtra      = string(T.extrainfo);

    % ---- Robust baseExtra parsing ----
    % Strip tokens like _R1, _R2, _R12, and _<number or decimal>kPa (case-insensitive)
    baseExtra = regexprep(rawExtra, '(?i)_?R\d+|_?\d+(\.\d+)?kPa', '');
    baseExtra = regexprep(baseExtra, '_+', '_');
    baseExtra = strip(baseExtra, '_');

    % ---- Modifiers flag (used to exclude from baseline) ----
    % Keep the simple contains test for familiarity, but make it case-insensitive and generic for kPa:
    isModifier = false(size(rawExtra));
    isModifier = isModifier | contains(rawExtra, {'R1','R2'}, 'IgnoreCase', true);
    isModifier = isModifier | ~cellfun(@isempty, regexpi(rawExtra, '\d+(\.\d+)?kPa')); % any kPa-like token

    % ---- CT flag (tolerant) ----
    dc = upper(strtrim(string(T.diseaseCat)));
    isCT = ismember(dc, ["CT","CONTROL"]);   % extend if needed

    % ================= Baseline computation =================
    % Group always includes DataFolder + sampletypeCat + springCat
    if combinePureCTAcrossBaseExtra
        % Combine pure CT across baseExtra
        [grpID, gDataFolder, gSample, gSpring] = findgroups( ...
            dataFolderStr, T.sampletypeCat, T.springCat);

        % Compute pure CT mean per (DataFolder, sampletype, spring)
        ctMeanPerGroup = splitapply(@(vals, ctFlag, modFlag) ...
            mean(vals(ctFlag & ~modFlag), 'omitnan'), ...
            T.(metricField), isCT, isModifier, grpID);

        % Also compute how many pure CT rows contributed (for diagnostics)
        ctCountPerGroup = splitapply(@(ctFlag, modFlag) ...
            sum(ctFlag & ~modFlag), isCT, isModifier, grpID);

        % Map baseline to rows
        T.(ctMeanVar) = ctMeanPerGroup(grpID);

        % Normalize
        T.(normVar) = T.(metricField) ./ T.(ctMeanVar) * targetCT;

        % Summary
        BaselineSummary = table(gDataFolder, gSample, gSpring, ctCountPerGroup, ctMeanPerGroup, ...
            'VariableNames', {'DataFolder','sampletype','spring','nPureCT','ctMean'});

    else
        % Original behavior: separate baseline per baseExtra
        [grpID, gDataFolder, gSample, gSpring, gBaseExtra] = findgroups( ...
            dataFolderStr, T.sampletypeCat, T.springCat, baseExtra);

        ctMeanPerGroup = splitapply(@(vals, ctFlag, modFlag) ...
            mean(vals(ctFlag & ~modFlag), 'omitnan'), ...
            T.(metricField), isCT, isModifier, grpID);

        ctCountPerGroup = splitapply(@(ctFlag, modFlag) ...
            sum(ctFlag & ~modFlag), isCT, isModifier, grpID);

        T.(ctMeanVar) = ctMeanPerGroup(grpID);
        T.(normVar)   = T.(metricField) ./ T.(ctMeanVar) * targetCT;

        BaselineSummary = table(gDataFolder, gSample, gSpring, gBaseExtra, ctCountPerGroup, ctMeanPerGroup, ...
            'VariableNames', {'DataFolder','sampletype','spring','baseExtra','nPureCT','ctMean'});
    end

    % ---- Diagnostics ----
    missingBase = isnan(T.(ctMeanVar));
    if any(missingBase)
        warning('Normalization:%s', ...
            sprintf('%d rows lack a pure CT baseline in their group (%s set to NaN).', ...
            sum(missingBase), normVar));
    end

    % Pretty print summary
    disp('--- CT baselines ---');
    disp(BaselineSummary);
end


% function [T, BaselineSummary] = NormalizeToCTBaseline(T, metricField, targetCT)
% %NORMALIZETOCTBASELINE Normalize a metric to per-DataFolder CT baseline (join-free).
% %   [T, BaselineSummary] = NormalizeToCTBaseline(T)
% %   [T, BaselineSummary] = NormalizeToCTBaseline(T, metricField, targetCT)
% %
% % Inputs:
% %   T           - table with variables:
% %                 DataFolder, extrainfo, diseaseCat, sampletypeCat, springCat,
% %                 and the metric to normalize (default: 'FRETav').
% %   metricField - (char/string) name of metric field in T to normalize (default: 'FRETav').
% %   targetCT    - scalar value to scale CT mean to (default: 1; e.g., 1 or 100).
% %
% % Outputs:
% %   T               - table with added columns:
% %                     T.(<metricField>)ctMeanNorm  (per-group pure CT mean)
% %                     T.(<metricField>)Norm    (normalized metric scaled to targetCT)
% %   BaselineSummary - compact summary of baselines used, per group, with
% %                     its baseline column named (<metricField>)ctMean.
% %
% % Mirrors the original script:
% %   - Uses DataFolder (string), detects modifiers in 'extrainfo'.
% %   - Groups by (DataFolder, sampletypeCat, springCat, baseExtra).
% %   - Computes pure-CT mean per group, excluding rows with modifiers.
% %   - Normalizes all rows: metric / ctMean * targetCT.
% %   - Warns if any rows lack a CT baseline (sets FRETnormValue to NaN for those).
% 
%     % ---- Defaults ----
%     if nargin < 2 || isempty(metricField), metricField = 'FRETav'; end
%     if nargin < 3 || isempty(targetCT),    targetCT    = 1;        end
% 
%     % Normalize metricField to string
%     metricField = string(metricField);
% 
%     % ---- Dynamic variable names (validated) ----
%     ctMeanVar = matlab.lang.makeValidName(metricField + "ctMeanNorm");
%     normVar   = matlab.lang.makeValidName(metricField + "Norm");
% 
%     % ---- Basic validation (lightweight, but helpful) ----
%     mustHave = {'DataFolder','extrainfo','diseaseCat','sampletypeCat','springCat', metricField};
%     for i = 1:numel(mustHave)
%         if ~ismember(mustHave{i}, T.Properties.VariableNames)
%             error('NormalizeToCTBaseline:MissingVar', ...
%                   'Variable "%s" not found in T.', mustHave{i});
%         end
%     end
%     if ~isnumeric(T.(metricField))
%         error('NormalizeToCTBaseline:ValueNotNumeric', ...
%               'The metric "%s" must be numeric.', metricField);
%     end
%     if ~iscategorical(T.diseaseCat)
%         warning('NormalizeToCTBaseline:Type', 'T.diseaseCat is expected categorical; converting.');
%         T.diseaseCat = categorical(T.diseaseCat);
%     end
%     if ~iscategorical(T.sampletypeCat)
%         warning('NormalizeToCTBaseline:Type', 'T.sampletypeCat is expected categorical; converting.');
%         T.sampletypeCat = categorical(T.sampletypeCat);
%     end
%     if ~iscategorical(T.springCat)
%         warning('NormalizeToCTBaseline:Type', 'T.springCat is expected categorical; converting.');
%         T.springCat = categorical(T.springCat);
%     end
% 
%     % ================= Normalization to per-DataFolder CT baseline =================
%     % Use DataFolder instead of date to group together
%     dataFolderStr = string(T.DataFolder);
% 
%     % Grab extrainfo and detect modifiers
%     rawExtra    = string(T.extrainfo);       % e.g., 'si48_4h', 'si72', 'si72_1kPa'
%     modifierTags = {'R1','R2','1kPa','10kPa'};
%     isModifier   = contains(rawExtra, modifierTags);
% 
%     % Build "baseExtra" by removing modifiers (so 'si72_1kPa' → 'si72'; 'si48_4h' stays 'si48_4h')
%     baseExtra = regexprep(rawExtra, '(_?(R1|R2|1kPa|10kPa))', '');
%     baseExtra = strrep(baseExtra, '__', '_');
%     baseExtra = strip(baseExtra, '_');
% 
%     % Flags
%     isCT  = (T.diseaseCat == 'CT');
%     % isCCM = (T.diseaseCat == 'CCM');  % not required for normalization
% 
%     % Group definition:
%     % Same DataFolder + same sample type (TS/TL) + same spring + same baseExtra (si/time without modifiers)
%     [grpID, grpLevels_dataFolder, grpLevels_sampletype, grpLevels_spring, grpLevels_baseExtra] = ...
%         findgroups(dataFolderStr, T.sampletypeCat, T.springCat, baseExtra);
% 
%     % Compute PURE-CT baseline per group (exclude modifiers)
%     ctMeanPerGroup = splitapply(@(vals, ctFlag, modFlag) ...
%         mean(vals(ctFlag & ~modFlag), 'omitnan'), ...
%         T.(metricField), isCT, isModifier, grpID);
% 
%     % Map per-group baseline back to each row using grpID (no join needed)
%     T.(ctMeanVar) = ctMeanPerGroup(grpID);
% 
%     % Normalize all rows to their group's pure-CT mean
%     T.(normVar) = T.(metricField) ./ T.(ctMeanVar) * targetCT;
% 
%     % Diagnostics: flag rows without a baseline
%     missingBase = isnan(T.(ctMeanVar));
%     if any(missingBase)
%         warning('Normalization:%s', ...
%             sprintf('%d rows lack a pure CT baseline in their DataFolder (%s set to NaN).', ...
%             sum(missingBase), normVar));
%     end
% 
%     % (Optional) Compact summary of baselines used
%     BaselineSummary = table( ...
%         grpLevels_dataFolder, grpLevels_sampletype, grpLevels_spring, grpLevels_baseExtra, ctMeanPerGroup, ...
%         'VariableNames', {'DataFolder','sampletype','spring','baseExtra', char(ctMeanVar)} );
% 
%     % Display summary (kept from your script)
%     disp('--- CT baselines per (DataFolder, sampletype, spring, baseExtra) ---');
%     disp(BaselineSummary);
% end
