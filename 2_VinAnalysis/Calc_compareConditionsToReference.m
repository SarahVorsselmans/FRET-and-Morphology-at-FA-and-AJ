function comparisonTable = Calc_compareConditionsToReference(inputTable, metrics, options)
% CALC_COMPARECONDITIONSTOREFERENCE
%   Compares all unique conditions in inputTable to a reference condition
%   across chosen metrics, using aggregation + multiple deviation measures.
%
% INPUTS:
%   inputTable  - one table from CatTables (e.g. CatTables.CatAll)
%   metrics     - cell array of metric column names, e.g. {'AverageCount','AverageSize'}
%   refDisease  - reference diseaseCat value, e.g. 'CT'
%   refExtra    - reference extrainfoCat value, e.g. 'rest'
%
% OUTPUT:
%   comparisonTable - table with one row per non-reference condition,
%                     columns: diseaseCat, extrainfoCat, conditionLabel,
%                     then per metric: _mean, _mean_ref, _absDiff,
%                     _relDiff_pct, _log2FC, _cohensD

arguments
    inputTable              table
    metrics                 cell    = {}
    options.refDisease      string = 'CT'
    options.refExtra        string = 'rest'
end
refDisease = categorical(options.refDisease);
refExtra   = categorical(options.refExtra);

    % --- 1. Find all unique conditions ---
    allCombos = unique(inputTable(:, {'diseaseCat', 'extrainfoCat'}));

    % --- 2. Isolate reference rows ---
    inputTable.diseaseCat = categorical(inputTable.diseaseCat);     % sometimes they are categorical & ordinal
    inputTable.extrainfoCat = categorical(inputTable.extrainfoCat); % sometimes they are categorical & ordinal
    isRef = inputTable.diseaseCat == refDisease & ...
            inputTable.extrainfoCat == refExtra;

    refRows = inputTable(isRef, :);

    if isempty(refRows)
        error('Reference condition (diseaseCat="%s", extrainfoCat="%s") not found.', ...
              refDisease, refExtra);
    end

    % --- 3. Pre-compute reference stats per metric ---
    refMean = zeros(1, numel(metrics));
    refStd  = zeros(1, numel(metrics));
    refN    = height(refRows);
    for m = 1:numel(metrics)
        vals = refRows.(metrics{m});
        refMean(m) = mean(vals, 'omitnan');
        refStd(m)  = std(vals,  'omitnan');
    end

    % --- 4. Loop over non-reference conditions ---
    % Remove reference from combo list
    allCombos.diseaseCat = categorical(allCombos.diseaseCat);     % sometimes they are categorical & ordinal
    allCombos.extrainfoCat = categorical(allCombos.extrainfoCat); % sometimes they are categorical & ordinal
    isRefCombo = allCombos.diseaseCat == refDisease & ...
                 allCombos.extrainfoCat == refExtra;
    condCombos = allCombos(~isRefCombo, :);
    % Check
    % disp(allCombos)
    % disp(isRefCombo)

    nConds   = height(condCombos);
    nMetrics = numel(metrics);

    % Pre-allocate output arrays
    condDiseaseCol = cell(nConds, 1);
    condExtraCol   = cell(nConds, 1);
    condLabelCol   = cell(nConds, 1);

    % One matrix per deviation type  [nConds x nMetrics]
    meanCond   = NaN(nConds, nMetrics);
    absDiff    = NaN(nConds, nMetrics);
    relDiff    = NaN(nConds, nMetrics);
    log2FC     = NaN(nConds, nMetrics);
    cohensD    = NaN(nConds, nMetrics);

    for c = 1:nConds
        cd = char(condCombos.diseaseCat(c));
        ce = char(condCombos.extrainfoCat(c));

        condDiseaseCol{c} = cd;
        condExtraCol{c}   = ce;
        condLabelCol{c}   = sprintf('%s_%s', cd, ce);

        % Rows for this condition
        isCond = inputTable.diseaseCat == condCombos.diseaseCat(c) & ...
                 inputTable.extrainfoCat == condCombos.extrainfoCat(c);
        condRows  = inputTable(isCond, :);
        condN     = height(condRows);

        for m = 1:nMetrics
            vals    = condRows.(metrics{m});
            cMean   = mean(vals, 'omitnan');
            cStd    = std(vals,  'omitnan');

            meanCond(c, m) = cMean;
            absDiff(c, m)  = cMean - refMean(m);

            % Relative difference (guard against zero reference)
            if refMean(m) ~= 0
                relDiff(c, m) = (cMean - refMean(m)) / abs(refMean(m)) * 100;
                log2FC(c, m)  = log2(cMean / refMean(m));
            end

            % Cohen's d with pooled SD
            pooledSD = sqrt(((condN - 1) * cStd^2 + (refN - 1) * refStd(m)^2) ...
                            / (condN + refN - 2));
            if pooledSD > 0
                cohensD(c, m) = (cMean - refMean(m)) / pooledSD;
            end
        end
    end

    % --- 5. Assemble output table ---
    comparisonTable = table(condDiseaseCol, condExtraCol, condLabelCol, ...
        'VariableNames', {'diseaseCat', 'extrainfoCat', 'conditionLabel'});

    for m = 1:nMetrics
        mn = metrics{m};  % metric name shorthand
        comparisonTable.(sprintf('%s_mean',       mn)) = meanCond(:, m);
        comparisonTable.(sprintf('%s_mean_ref',   mn)) = repmat(refMean(m), nConds, 1);
        comparisonTable.(sprintf('%s_absDiff',    mn)) = absDiff(:,  m);
        comparisonTable.(sprintf('%s_relDiff_pct',mn)) = relDiff(:,  m);
        comparisonTable.(sprintf('%s_log2FC',     mn)) = log2FC(:,   m);
        comparisonTable.(sprintf('%s_cohensD',    mn)) = cohensD(:,  m);
    end
end