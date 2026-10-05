function hFig = PlotFRETTSCCMBoxplot(T, lineval, valueVar, includePatterns, excludePatterns) %, includePatterns, excludePatterns)
%PLOTFRETBOXPLOT Minimal boxplot visualization mirrored from your script.
%   hFig = PlotFRETBoxplot(T)
%   hFig = PlotFRETBoxplot(T, valueVar)
%
% Inputs:
%   T        - table with variables: xCat (categorical), disease (cellstr/strings),
%              sampletypeCat (categorical or char), date (string/categorical),
%              and the numeric value variable (e.g., 'FRETav' or 'FRETavNorm').
%   valueVar - name of the value variable in T to plot (default: 'FRETav').
%
% Output:
%   hFig     - figure handle.

% Ordere of boxplots is determined in the main code under "Build sorting table"

    % if nargin < 2 || isempty(valueVar)
    %     valueVar = 'FRETav'; % default to your current field
    % end
    % if nargin < 3
    %     Info = struct(); % allow calling without Info, but ordering will be whatever T.xCat has
    % end

    % % ----------------------- Filter data input ----------------------- % this doesn't works because categories stay remembered
    % % includePatterns = {'VinTS_FL'}; % careful: include with ... or with ...
    % % excludePatterns = {'ish'}; %'R1', 'R2', 'kPa'
    % includeMask = false(height(T),1);
    % for i = 1:numel(includePatterns)
    %     includeMask = includeMask | contains(T.SampleFolder, includePatterns{i});
    % end
    % 
    % excludeMask = false(height(T),1);
    % for i = 1:numel(excludePatterns)
    %     excludeMask = excludeMask | contains(T.SampleFolder, excludePatterns{i});
    % end
    % 
    % T = T(includeMask & ~excludeMask, :);

    % Filter
    T = Filter_IncludeExclude(T, includePatterns, excludePatterns);

    % ---- Basic validation (minimal, but helpful) ----
    mustHave = {'xCat', 'disease', 'sampletypeCat', 'date', valueVar};
    for i = 1:numel(mustHave)
        if ~ismember(mustHave{i}, T.Properties.VariableNames)
            error('PlotFRETBoxplot:MissingVar', ...
                  'Variable "%s" not found in T.', mustHave{i});
        end
    end
    if ~iscategorical(T.xCat)
        error('PlotFRETBoxplot:Type', 'T.xCat must be categorical.');
    end
    if ~isnumeric(T.(valueVar))
        error('PlotFRETBoxplot:Type', 'T.%s must be numeric.', valueVar);
    end

    % ===== DATE-BASED ORDERING using T.DataFolder =====
    % Convert to string to make sorting safe
    df = string(T.DataFolder);
    % Sort table by DataFolder so that identical folders stay grouped
    [T, ~] = sortrows(T, 'DataFolder');
    % After sorting, rebuild xCat in the same order
    T.xCat = removecats(categorical(T.xCat, 'Ordinal', true));
    % Extract the grouping categories in their actual order of appearance
    groupCats = cellstr(unique(T.xCat, 'stable'));

    % ================= MAIN: Boxplot visualization =================
    % ---------- Create boxplot ----------
    hFig = figure('Color','w'); %hold on;
    ax = axes('Parent', hFig, 'Color', 'w'); % w: ensure white axes background
    hold(ax, 'on');

    % Boxplot using selected value variable vs grouping xCat
    %groupCats = categories(removecats(T.xCat)); % get only used categories
    boxWidth = 0.6;   % for jitter?
    boxplot(T.(valueVar), T.xCat, ...
        'GroupOrder', groupCats, ...
        'Whisker', 1.5, ...
        'Symbol', 'k.', ...
        'Widths', boxWidth, ...
        'Colors', 'k', ...
        'Parent',  ax);  % k: base color; we’ll recolor patches/medians below

    set(ax,'XTickLabelRotation',45);

    % Add jitter
    % Prepare jittered x-positions that align with box centers
    % Enforce the same category order before mapping to numeric indices
    T.xCatJit = categorical(T.xCat, groupCats, 'Ordinal', true);
    T.xCatJit = reordercats(T.xCatJit, groupCats);
    xIdx = grp2idx(T.xCatJit);          % 1..K for each row
    yVal = T.(valueVar);
    % Optional: drop NaNs so scatter sizes match
    mask = ~isnan(yVal);
    xIdx = xIdx(mask);
    yVal = yVal(mask);
    % Jitter amplitude (keep inside the box width)
    % ~80% of half-width looks good visually
    jitterAmp = (boxWidth/2) * 0.8;  
    xJit = xIdx + (rand(size(xIdx)) - 0.5) * 2 * jitterAmp;
    % Overlay jittered scatter
    scatter(ax, xJit, yVal, ...
        18, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.5, ...
        'MarkerEdgeColor', 'none');

    % Add horizontal line at y = 28
    yline(ax, lineval, '--g', 'LineWidth', 1.5);

    % ---------- Style boxes ----------
    hBox = findobj(ax, 'Type','line', 'Tag','Box');
    xData = get(hBox,'XData'); if ~iscell(xData), xData = {xData}; end
    xCenter = cellfun(@(x) mean(x(:)), xData);      % mean x of top edge
    [~, idxL2R_B] = sort(xCenter, 'ascend');
    hBox = hBox(idxL2R_B);
    
    % Use visible tick labels as group keys (matches left→right)
    %groupNames = categories(T.xCat);

    for i = 1:numel(hBox)
        thisGroup = groupCats{i};

        % representative row for this SampleFolder
        rowIdx = find(T.xCat == thisGroup, 1, 'first');

        % Disease color
        if iscell(T.disease)
            thisDisease = string(T.disease{rowIdx});
        else
            thisDisease = string(T.disease(rowIdx));
        end
        if strcmp(thisDisease,'CT')
            col = [0 0 1];
        else
            col = [1 0 0];
        end

        % Sample type: TS filled, TL outline
        stype = T.sampletypeCat(rowIdx);
        if ~(ischar(stype) || isstring(stype))
            stype = string(stype);
        end
        isTS = strcmp(stype, 'TS');

        % Redraw box
        x = get(hBox(i),'XData');
        y = get(hBox(i),'YData');
        set(hBox(i),'Visible','off'); % hide original
        if isTS
            patch('XData', x, 'YData', y, 'FaceColor', col, 'FaceAlpha', 0.6, ...
                  'EdgeColor', col, 'LineWidth', 1.5, 'Parent', ax);
        else
            patch('XData', x, 'YData', y, 'FaceColor', 'none', ...
                  'EdgeColor', col, 'LineWidth', 1.5, 'Parent', ax);

        % % Access disease and type
        % if iscell(T.disease)
        %     thisDisease = T.disease{rowIdx};
        % else
        %     thisDisease = string(T.disease(rowIdx));
        % end
        % thisType = T.sampletypeCat(rowIdx);
        % 
        % % disease color
        % if strcmp(thisDisease,'CT')
        %     col = [0 0 1];   % blue
        % else
        %     col = [1 0 0];   % red
        % end
        % 
        % % get box coordinates
        % x = get(hBox(i),'XData');
        % y = get(hBox(i),'YData');
        % 
        % % remove original box
        % set(hBox(i),'Visible','off');
        % 
        % % filled for TS, outline for TL
        % if ischar(thisType) || isstring(thisType)
        %     isTS = strcmp(thisType, 'TS');
        % else
        %     % categorical comparison
        %     isTS = (thisType == 'TS');
        % end
        % 
        % if isTS
        %     % filled box
        %     patch(x, y, col, ...
        %         'FaceAlpha',0.6, ...
        %         'EdgeColor',col, ...
        %         'LineWidth',1.5);
        % else
        %     % TL: outline only
        %     patch(x, y, col, ...
        %         'FaceColor','none', ...
        %         'EdgeColor',col, ...
        %         'LineWidth',1.5);
        end
    end

    % ---------- Recolor median lines to match disease ----------
    hMedian = findobj(ax, 'Type','line', 'Tag','Median');
    xm = get(hMedian,'XData'); if ~iscell(xm), xm = {xm}; end
    xc = cellfun(@(x) mean(x(:)), xm);
    [~, idxL2R_M] = sort(xc,'ascend');
    hMedian = hMedian(idxL2R_M);   % left-to-right, same as boxes

    % groupNames = categories(T.xCat);
    for i = 1:numel(hMedian)
        grp = groupCats{i};
        rowIdx = find(T.xCat == grp, 1, 'first');

        if iscell(T.disease)
            thisDisease = string(T.disease{rowIdx});
        else                
            thisDisease = string(T.disease(rowIdx));
        end
        if strcmp(thisDisease,'CT')
            col = [0 0 1];
        else
            col = [1 0 0];
        end

        set(hMedian(i), 'Color', col, 'LineWidth', 1.8);

    % for i = 1:numel(hMedian)
    %     thisGroup = groupNames{i};
    % 
    %     % representative row for this SampleFolder
    %     rowIdx = find(T.xCat == thisGroup, 1, 'first');
    % 
    %     if iscell(T.disease)
    %         thisDisease = T.disease{rowIdx};
    %     else
    %         thisDisease = string(T.disease(rowIdx));
    %     end
    % 
    %     % disease color
    %     if strcmp(thisDisease,'CT')
    %         col = [0 0 1];   % blue
    %     else
    %         col = [1 0 0];   % red
    %     end
    % 
    %     set(hMedian(i), ...
    %         'Color', col, ...
    %         'LineWidth', 1.8);
    end

    % ---------- Add vertical dotted lines between date groups ----------
    %dates = T.date;
    dfPerBox = strings(numel(groupCats),1);
    for i = 1:numel(groupCats)
        rowIdx = find(T.xCat == groupCats{i}, 1, 'first');
        %datePerBox(i) = string(T.date(rowIdx));
        dfPerBox(i) = string(T.DataFolder(rowIdx));
    end
    % Find boundaries where DataFolder changes
    %changeIdx = find(~strcmp(datePerBox(1:end-1), datePerBox(2:end)));
    changeIdx = find(~strcmp(dfPerBox(1:end-1), dfPerBox(2:end)));
    % Draw vertical dotted lines between boxes
    for k = changeIdx(:).'
        xline(ax, k+0.5, ':k', 'LineWidth', 1);
    end

    % ---------- Axis tick formatting ----------
    %ax = gca;

    % Remove top and right axes + ticks
    set(ax, ...
        'Box','off', ...
        'TickDir','out', ...
        'XAxisLocation','bottom', ...
        'YAxisLocation','left', ...
        'XMinorTick','off', ...
        'YMinorTick','off');

    % ---------- Legend (manual) ----------
    %hold on;
    p1 = plot(ax, nan,nan,'s','MarkerFaceColor','b','MarkerEdgeColor','b'); % CT
    p2 = plot(ax, nan,nan,'s','MarkerFaceColor','r','MarkerEdgeColor','r'); % CCM
    p3 = plot(ax, nan,nan,'s','MarkerEdgeColor','k','MarkerFaceColor','none'); % TL

    legend(ax, [p1 p2 p3], ...
        {'CT','CCM','TL'}, ...
        'Location','northeast');

    % Y-axis settings: keep your defaults; adjust if valueVar suggests normalized
    if strcmpi(valueVar, 'FRETavNorm')
        ylim(ax, [0.8 1.2]);
        ylabel(ax, 'Normalized FRET value');
    elseif strcmpi(valueVar, 'FRETav')
        ylim(ax, [20 40]);
        ylabel(ax, 'FRET efficiency (%)');
    else
        %ylim(ax, 'auto');
        ylabel(ax, valueVar, 'Interpreter','none');
    end

    title(ax, 'FRET-TS efficiency, grouped by date');

    hold(ax, 'off');
end
