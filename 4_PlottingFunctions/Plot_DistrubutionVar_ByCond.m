function Plot_DistrubutionVar_ByCond(T, innerVar, options)
% Plot_DistrubutionVar_ByCond
% --------------------------------------------------------------
% Plots the distribution of a variable by condition (diseaseCat x extrainfoCat).
%
% INPUTS:
%   T         : table containing FA metrics
%               - either as cell arrays of numeric vectors per row, OR
%               - as a simple numeric value per row (scalar)
%   innerVar  : name of column to plot (e.g., "Ecc", "MajAx")
%   nBins     : histogram bins (default: 40)
%   doClip    : clip 1–99 percentile for plotting (default: true)
%
%   Optional filtering based on another FA metric:
%   filterVar : name of variable to filter on (same shape logic as innerVar)
%   filterOp  : '<' '<=' '>' '>='
%   filterVal : numeric threshold
%
%   levels    : (optional) cell array of strings — ordered subset of
%               extrainfoCat values to include, e.g. {'rest','R1','R2'}.
%               If empty or omitted, all unique values found in T.extrainfoCat
%               are used in order of first appearance.
%
% Groups:
%   T.diseaseCat   ∈ {CT, CCM}
%   T.extrainfoCat ∈ dynamically resolved from T (or restricted to `levels`)
% --------------------------------------------------------------

arguments
    T                               table
    innerVar                        (1,1) string
    options.nBins                   (1,1) double  = 40 % if empty -> value 40
    options.doClip                  (1,1) logical = true
    options.levels                  cell           = {}

    options.filterVar1   (1,1) string = ""
    options.filterOp1    (1,1) string {mustBeMember(options.filterOp1, ["<","<=",">",">=",""])} = ""
    options.filterVal1   (1,1) double = NaN
    options.filterVar2   (1,1) string = ""
    options.filterOp2    (1,1) string {mustBeMember(options.filterOp2, ["<","<=",">",">=",""])} = ""
    options.filterVal2   (1,1) double = NaN
end

%% Rename to be easier to use
nBins     = options.nBins;
doClip    = options.doClip;
levels    = options.levels;

filterVar1 = options.filterVar1;
filterOp1  = options.filterOp1;
filterVal1 = options.filterVal1;
filterVar2 = options.filterVar2;
filterOp2  = options.filterOp2;
filterVal2 = options.filterVal2;
if filterVar1 == "", filterVar1 = []; end
if filterVar2 == "", filterVar2 = []; end

    %% -------- RESOLVE LEVELS FROM TABLE --------
    % Retrieve unique extrainfoCat values in order of first appearance
    if iscategorical(T.extrainfoCat)
        allLevels = cellstr(unique(T.extrainfoCat, 'stable'));
    else
        allLevels = unique(T.extrainfoCat, 'stable');
        if isstring(allLevels), allLevels = cellstr(allLevels); end
    end

    % If user supplied levels, validate and use that order; otherwise use all
    if ~isempty(levels)
        % Normalise to cell array of char
        if isstring(levels), levels = cellstr(levels); end
        badLevels = levels(~ismember(levels, allLevels));
        if ~isempty(badLevels)
            error('Requested level(s) not found in extrainfoCat: %s', ...
                  strjoin(badLevels, ', '));
        end
        infoList = levels;
    else
        infoList = allLevels;
    end
    nLevels = numel(infoList);

    diseaseList = ["CT","CCM"];

    %% -------- GLOBAL LIMITS (from ALL rows of T) --------
    allGlobalVals = [];

    for i = 1:height(T)
        vals = getRowVals(T.(innerVar), i);
        if isempty(vals), continue; end
        allGlobalVals = [allGlobalVals; vals(:)]; %#ok<AGROW>
    end

    allGlobalVals = allGlobalVals(~isnan(allGlobalVals));

    if isempty(allGlobalVals)
        warning('No valid numeric values found in column "%s". Nothing to plot.', innerVar);
        return;
    end

    if doClip
        xLow  = prctile(allGlobalVals, 1);
        xHigh = prctile(allGlobalVals, 99);
    else
        xLow  = min(allGlobalVals);
        xHigh = max(allGlobalVals);
    end

    %% -------- Plotting --------
    figure;
    tiledlayout(2, nLevels, 'TileSpacing', 'compact', 'Padding', 'compact');

    for d = 1:numel(diseaseList)
        for e = 1:nLevels

            nexttile;

            mask = T.diseaseCat == diseaseList(d) & ...
                   T.extrainfoCat == infoList{e};

            subT = T(mask, :);

            allVals = [];

            for i = 1:height(subT)
                plotVals = getRowVals(subT.(innerVar), i);
                if isempty(plotVals)
                    continue
                end
            
                % ----------------------------------------------------
                % Build one combined mask from all filters
                % ----------------------------------------------------
            
                if isscalar(plotVals)
                    keep = true;
                else
                    keep = true(size(plotVals));
                end
            
                %% ----- Filter 1 -----
            
                if ~isempty(filterVar1)
                    fVals = getRowVals(subT.(filterVar1), i);
                    
                    if isempty(fVals)
                        continue
                    end
            
                    if isscalar(fVals)
                        keep = keep & ...
                            applyOpScalar(fVals, filterVal1, filterOp1);
            
                    else
                        if isscalar(plotVals)
                            keep = keep & ...
                                any(applyOpElementwise(fVals, filterVal1, filterOp1));
            
                        else
                            if numel(fVals) ~= numel(plotVals)
                                warning(['Row %d: filterVar "%s" and innerVar "%s" ' ...
                                         'lengths differ (%d vs %d). Skipping row.'], ...
                                         i, filterVar1, innerVar, ...
                                         numel(fVals), numel(plotVals));
                                continue
                            end
                            keep = keep & ...
                                applyOpElementwise(fVals, filterVal1, filterOp1);
                        end
                    end
                end
            
                %% ----- Filter 2 -----
            
                if ~isempty(filterVar2)
                    fVals = getRowVals(subT.(filterVar2), i);
            
                    if isempty(fVals)
                        continue
                    end
            
                    if isscalar(fVals)
                        keep = keep & ...
                            applyOpScalar(fVals, filterVal2, filterOp2);
                    else
            
                        if isscalar(plotVals)
                            keep = keep & ...
                                any(applyOpElementwise(fVals, filterVal2, filterOp2));
                        else
                            if numel(fVals) ~= numel(plotVals)
            
                                warning(['Row %d: filterVar "%s" and innerVar "%s" ' ...
                                         'lengths differ (%d vs %d). Skipping row.'], ...
                                         i, filterVar2, innerVar, ...
                                         numel(fVals), numel(plotVals));
                                continue
                            end
                            keep = keep & ...
                                applyOpElementwise(fVals, filterVal2, filterOp2);
                        end
                    end
                end
            
                %% ----- Apply combined mask once -----
                if isscalar(plotVals)
                    if ~keep
                        continue
                    end
            
                else
                    plotVals = plotVals(keep);
                    if isempty(plotVals)
                        continue
                    end
                end
            
                allVals = [allVals; plotVals(:)]; %#ok<AGROW>
            end

            if doClip
                allVals = allVals(allVals >= xLow & allVals <= xHigh);
            end

            if ~isempty(allVals)
                histogram(allVals, nBins);
                xlim([xLow xHigh]);

                meanVal = mean(allVals, 'omitnan');
                n       = numel(allVals);
                if ~isnan(meanVal)
                    text(0.98, 0.95, sprintf('Av = %.3g \n N = %d', meanVal, n), ...
                        'Units',               'normalized', ...
                        'HorizontalAlignment', 'right', ...
                        'VerticalAlignment',   'top', ...
                        'FontWeight',          'bold', ...
                        'Color',               [0.2 0.2 0.2]);
                end
            else
                text(0.5, 0.5, 'No data', ...
                    'HorizontalAlignment', 'center', ...
                    'Units',               'normalized', ...
                    'VerticalAlignment',   'middle');
                xlim([0 1]); ylim([0 1]);
            end

            title(sprintf('%s - %s', diseaseList(d), infoList{e}));
            xlabel(string(innerVar), 'Interpreter', 'none');
            ylabel('Abundance');
            grid on;

        end
    end

    %% ---- Figure-level title ----
    titleTxt = sprintf('Distribution of %s', string(innerVar));
    if ~isempty(filterVar1)
        titleTxt = sprintf('%s | %s %s %.4g', ...
            titleTxt, filterVar1, filterOp1, filterVal1);
    end
    if ~isempty(filterVar2)
        titleTxt = sprintf('%s | %s %s %.4g', ...
            titleTxt, filterVar2, filterOp2, filterVal2);
    end
    sgtitle(titleTxt, 'Interpreter', 'none');
end

% ---------- Helper: get row values from a column that can be cell-of-vectors or numeric ----------
function vals = getRowVals(col, idx)
    if iscell(col)
        vals = col{idx};
        if isempty(vals), return; end
        if ~isnumeric(vals) || (~isvector(vals) && ~isscalar(vals))
            error('Row %d contains non-numeric or non-vector/scalar data.', idx);
        end
    elseif isnumeric(col)
        if size(col, 2) ~= 1
            error('Numeric column must be Nx1.');
        end
        vals = col(idx, 1);
    else
        error('Column type must be cell array (of numeric vectors/scalars) or numeric Nx1.');
    end
end

% ---------- Helper: filters ----------
function plotVals = applyFilterToRow( ...
            plotVals, fVals, filterVal, filterOp, ...
            rowIdx, filterVar, innerVar)
    if isempty(fVals)
        plotVals = [];
        return
    end
    if isscalar(fVals) && ~isscalar(plotVals)
        keep = applyOpElementwise(fVals * ones(size(plotVals)), ...
                                  filterVal, filterOp);
        plotVals = plotVals(keep);
    elseif ~isscalar(fVals) && isscalar(plotVals)
        keepRow = any(applyOpElementwise(fVals, ...
                                        filterVal, filterOp));
        if ~keepRow
            plotVals = [];
        end
    elseif ~isscalar(fVals) && ~isscalar(plotVals)
        if numel(fVals) ~= numel(plotVals)
            warning(['Row %d: filterVar "%s" and innerVar "%s" lengths differ (%d vs %d). ' ...
                     'Skipping row.'], ...
                     rowIdx, filterVar, innerVar, ...
                     numel(fVals), numel(plotVals));
            plotVals = [];
        else
            keep = applyOpElementwise(fVals, ...
                                      filterVal, filterOp);
            plotVals = plotVals(keep);
        end
    else
        keepScalar = applyOpScalar(fVals, ...
                                   filterVal, filterOp);
        if ~keepScalar
            plotVals = [];
        end
    end
end

% ---------- Helper: elementwise comparison (vector vs threshold) ----------
function keep = applyOpElementwise(arr, thr, op)
    switch op
        case '<',   keep = arr <  thr;
        case '<=',  keep = arr <= thr;
        case '>',   keep = arr >  thr;
        case '>=',  keep = arr >= thr;
        otherwise,  error('Invalid filter operator.');
    end
end

% ---------- Helper: scalar comparison ----------
function keep = applyOpScalar(a, thr, op)
    switch op
        case '<',   keep = a <  thr;
        case '<=',  keep = a <= thr;
        case '>',   keep = a >  thr;
        case '>=',  keep = a >= thr;
        otherwise,  error('Invalid filter operator.');
    end
end

% ---------- Helper: broadcast scalar filter to a vector of plot values ----------
function keep = applyOpToArray(arr, fScalar, op, mode)
    switch mode
        case 'broadcastFilterScalar'
            error('Internal misuse: applyOpToArray requires explicit compare with threshold. Use applyOpElementwise/Scalar.');
        case 'compareToThreshold'
            error('Internal misuse: compareToThreshold should not be called with two args.');
        otherwise
            error('Unknown mode.');
    end
end

% function Plot_DistrubutionVar_ByCond(T, innerVar, nBins, doClip, filterVar, filterOp, filterVal)
% % Plot_DistrubutionVar_ByCond
% % --------------------------------------------------------------
% % Plots the distribution of a variable by condition (diseaseCat x extrainfoCat).
% %
% % INPUTS:
% %   T         : table containing FA metrics
% %               - either as cell arrays of numeric vectors per row, OR
% %               - as a simple numeric value per row (scalar)
% %   innerVar  : name of column to plot (e.g., "Ecc", "MajAx")
% %   nBins     : histogram bins (default: 40)
% %   doClip    : clip 1–99 percentile for plotting (default: true)
% %
% %   Optional filtering based on another FA metric:
% %   filterVar : name of variable to filter on (same shape logic as innerVar)
% %   filterOp  : '<' '<=' '>' '>='
% %   filterVal : numeric threshold
% %
% % Groups:
% %   T.diseaseCat   ∈ {CT, CCM}
% %   T.extrainfoCat ∈ {rest, R1, R2, 10kPa, 1kPa}
% % --------------------------------------------------------------
% 
%     if nargin < 3 || isempty(nBins),   nBins = 40; end
%     if nargin < 4 || isempty(doClip),  doClip = true; end
%     if nargin < 5, filterVar = []; end
% 
%     innerVar  = string(innerVar);
%     if ~isempty(filterVar), filterVar = string(filterVar); end
% 
%     % ---- Validate inputs ----
%     if ~istable(T)
%         error('First input must be a table.');
%     end
%     if ~ismember(innerVar, string(T.Properties.VariableNames))
%         error('Column "%s" not found in the table.', innerVar);
%     end
%     if ~all(ismember(["diseaseCat","extrainfoCat"], string(T.Properties.VariableNames)))
%         error('Table must contain columns "diseaseCat" and "extrainfoCat".');
%     end
%     if ~isempty(filterVar) && ~ismember(filterVar, string(T.Properties.VariableNames))
%         error('Filter column "%s" not found in the table.', filterVar);
%     end
%     if ~isempty(filterVar) && ~ismember(filterOp, ["<","<=",">",">="])
%         error('filterOp must be <, <=, >, or >=.');
%     end
% 
%     diseaseList = ["CT","CCM"];
%     infoList    = ["rest","R1","R2","10kPa","1kPa"];
% 
%     %% -------- GLOBAL LIMITS (from ALL rows of T) --------
%     allGlobalVals = [];
% 
%     for i = 1:height(T)
%         vals = getRowVals(T.(innerVar), i);
%         if isempty(vals), continue; end
%         allGlobalVals = [allGlobalVals; vals(:)]; %#ok<AGROW>
%     end
% 
%     allGlobalVals = allGlobalVals(~isnan(allGlobalVals));
% 
%     if isempty(allGlobalVals)
%         warning('No valid numeric values found in column "%s". Nothing to plot.', innerVar);
%         return;
%     end
% 
%     if doClip
%         xLow  = prctile(allGlobalVals,1);
%         xHigh = prctile(allGlobalVals,99);
%     else
%         xLow  = min(allGlobalVals);
%         xHigh = max(allGlobalVals);
%     end
% 
%     %% -------- Plotting --------
%     figure;
%     tiledlayout(2,5,'TileSpacing','compact','Padding','compact');
% 
%     for d = 1:numel(diseaseList)
%         for e = 1:numel(infoList)
% 
%             nexttile;
% 
%             mask = T.diseaseCat == diseaseList(d) & ...
%                    T.extrainfoCat == infoList(e);
% 
%             subT = T(mask,:);
% 
%             allVals = [];
% 
%             for i = 1:height(subT)
%                 plotVals = getRowVals(subT.(innerVar), i);
%                 if isempty(plotVals), continue; end
% 
%                 % ----- Optional filtering -----
%                 if ~isempty(filterVar)
%                     fVals = getRowVals(subT.(filterVar), i);
% 
%                     % Align shapes for filtering:
%                     % - If filter is scalar and plotVals is vector: broadcast scalar.
%                     % - If both vectors: lengths must match.
%                     % - If plotVals is scalar and filter is vector: if lengths differ, use row-level keep.
%                     if isempty(fVals)
%                         % No filter values for this row -> skip row entirely
%                         continue;
%                     end
% 
%                     if isscalar(fVals) && ~isscalar(plotVals)
%                         keep = applyOpToArray(plotVals, fVals, filterOp, 'broadcastFilterScalar');
%                     elseif ~isscalar(fVals) && isscalar(plotVals)
%                         % If filter vector but inner scalar, reduce vector filter to a single decision
%                         % (keep the scalar if ANY elements satisfy the condition).
%                         keepRow = any(applyOpToArray(fVals, filterVal, filterOp, 'compareToThreshold'));
%                         if ~keepRow
%                             plotVals = []; % drop row
%                         end
%                         % else keep the scalar as-is
%                     else
%                         % Both vectors or both scalars
%                         if ~isscalar(fVals) && ~isscalar(plotVals)
%                             if numel(fVals) ~= numel(plotVals)
%                                 warning(['Row %d: filterVar "%s" and innerVar "%s" lengths differ (%d vs %d). ' ...
%                                          'Skipping this row.'], i, filterVar, innerVar, numel(fVals), numel(plotVals));
%                                 plotVals = [];
%                             else
%                                 keep = applyOpElementwise(fVals, filterVal, filterOp);
%                                 plotVals = plotVals(keep);
%                             end
%                         else
%                             % both scalars
%                             keepScalar = applyOpScalar(fVals, filterVal, filterOp);
%                             if ~keepScalar
%                                 plotVals = [];
%                             end
%                         end
%                     end
% 
%                     % If we computed a 'keep' mask via broadcast case:
%                     if exist('keep','var') && ~isempty(plotVals) && ~isscalar(plotVals)
%                         plotVals = plotVals(keep);
%                         clear keep;
%                     end
%                 end
% 
%                 if isempty(plotVals), continue; end
% 
%                 allVals = [allVals; plotVals(:)];
%             end
% 
%             allVals = allVals(~isnan(allVals));
% 
%             if doClip
%                 allVals = allVals(allVals >= xLow & allVals <= xHigh);
%             end
% 
%             % if ~isempty(allVals)
%             %     histogram(allVals, nBins);
%             %     xlim([xLow xHigh]);
%             % else
%             %     text(0.5,0.5,'No data','HorizontalAlignment','center');
%             %     xlim([0 1]); ylim([0 1]); % ensure text is visible
%             % end
%             % 
%             if ~isempty(allVals)
%                 h = histogram(allVals, nBins);
%                 xlim([xLow xHigh]);
% 
%                 % Mean annotation
%                 meanVal = mean(allVals, 'omitnan');    % mean of values actually plotted
%                 n = numel(allVals);
%                 if ~isnan(meanVal)
%                     % Put label in the top-right corner of the axes, independent of data limits
%                     text(0.98, 0.95, sprintf('Av = %.3g \n N = %d', meanVal, n), ...
%                         'Units','normalized', ...
%                         'HorizontalAlignment','right', ...
%                         'VerticalAlignment','top', ...
%                         'FontWeight','bold', ...
%                         'Color',[0.2 0.2 0.2]);
%                 end
%             else
%                 text(0.5,0.5,'No data','HorizontalAlignment','center', ...
%                     'Units','normalized', 'VerticalAlignment','middle');
%                 xlim([0 1]); ylim([0 1]); % ensure text is visible
%             end
%             title(sprintf('%s - %s', diseaseList(d), infoList(e)));
%             xlabel(string(innerVar), 'Interpreter','none');
%             ylabel('Abundance');
%             grid on;
% 
%         end
%     end
% 
%     %% Title
%     if isempty(filterVar)
%         sgtitle(sprintf('Distribution of %s', string(innerVar)), 'Interpreter','none');
%     else
%         sgtitle(sprintf('%s where %s %s %.4g', ...
%             string(innerVar), string(filterVar), filterOp, filterVal), 'Interpreter','none');
%     end
% end
% 
% % ---------- Helper: get row values from a column that can be cell-of-vectors or numeric ----------
% function vals = getRowVals(col, idx)
%     if iscell(col)
%         vals = col{idx};
%         if isempty(vals), return; end
%         if ~isnumeric(vals) || (~isvector(vals) && ~isscalar(vals))
%             error('Row %d contains non-numeric or non-vector/scalar data.', idx);
%         end
%     elseif isnumeric(col)
%         % Expect Nx1 numeric column; take scalar at this row
%         if size(col,2) ~= 1
%             error('Numeric column must be Nx1.');
%         end
%         vals = col(idx,1);
%     else
%         error('Column type must be cell array (of numeric vectors/scalars) or numeric Nx1.');
%     end
% end
% 
% % ---------- Helper: elementwise comparison (vector vs threshold) ----------
% function keep = applyOpElementwise(arr, thr, op)
%     switch op
%         case '<',   keep = arr <  thr;
%         case '<=',  keep = arr <= thr;
%         case '>',   keep = arr >  thr;
%         case '>=',  keep = arr >= thr;
%         otherwise, error('Invalid filter operator.');
%     end
% end
% 
% % ---------- Helper: scalar comparison ----------
% function keep = applyOpScalar(a, thr, op)
%     switch op
%         case '<',   keep = a <  thr;
%         case '<=',  keep = a <= thr;
%         case '>',   keep = a >  thr;
%         case '>=',  keep = a >= thr;
%         otherwise, error('Invalid filter operator.');
%     end
% end
% 
% % ---------- Helper: broadcast scalar filter to a vector of plot values ----------
% % mode 'broadcastFilterScalar': fVals is scalar, compare it to filterVal to decide which plotVals to keep.
% % mode 'compareToThreshold': compare array to threshold (used when reducing filter vector)
% function keep = applyOpToArray(arr, fScalar, op, mode)
%     switch mode
%         case 'broadcastFilterScalar'
%             % Decide keep based on comparing the scalar filter value to filterVal, then keep all or none
%             % But more useful: apply threshold using the scalar filter -> we need filterVal in scope.
%             % Since we don't have it here, redirect: keep all if scalar passes; else keep none.
%             error('Internal misuse: applyOpToArray requires explicit compare with threshold. Use applyOpElementwise/Scalar.');
%         case 'compareToThreshold'
%             error('Internal misuse: compareToThreshold should not be called with two args.');
%         otherwise
%             error('Unknown mode.');
%     end
% end
