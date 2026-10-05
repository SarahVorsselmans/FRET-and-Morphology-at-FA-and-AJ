function fig = Plot_FAClassDistribution(T, options)
% Plot_FAClassDistribution  Violin + scatter of FRET or Intensity per FA
% class, split by diseaseCat (CT/CCM) × extrainfoCat (rest/R1/R2/1kPa/10kPa).
% Scatter dot shape varies by DataFolder (experiment).
% Requires Calc_FAIntClassComparison to have been run first.
%
% USAGE:
%   fig = Plot_FAClassDistribution(T);
%   fig = Plot_FAClassDistribution(T, Metric='Int');
%   fig = Plot_FAClassDistribution(T, Metric='Int', XLimit=500);
%   fig = Plot_FAClassDistribution(T, ShowStats=false);
%
% DEPENDENCIES (must be on MATLAB path)
%   Filter_IncludeExclude.m
%   getPlotMarkers.m

arguments
    T          table
    options.Metric      (1,1) string {mustBeMember(options.Metric, ...
                                 ["FRET","Int","Frac"])} = "FRET"
    options.XLimit      (1,1) double  = Inf
    options.Levels      cell          = {}
    options.includePatterns cell      = {}
    options.excludePatterns cell      = {}
    options.legendNames cell          = {}
    options.ShowStats   (1,1) logical = true   % overlay ANOVA + sig brackets
    options.Alpha       (1,1) double  = 0.05   % significance threshold
end

% extract number of categories
%nCats = numel(categories);

    %% Filter rows via include/exclude on SampleFolder
    T = Filter_IncludeExclude(T, options.includePatterns, options.excludePatterns);
    
    %% ── 1. CHECK PRE-COMPUTED COLUMNS ────────────────────────────────────────
    if options.Metric == "FRET"
        required = {'FA_summaryFRETPx'};
    elseif options.Metric == "Int"
        required = {'FA_summaryIntPx'};
    else
        required = {'FA_summaryFrac'};
    end
    missing = required(~ismember(required, T.Properties.VariableNames));
    if ~isempty(missing)
        error('Missing columns: %s\nRun Calc_FAIntClassComparison(T) first.', ...
              strjoin(missing, ', '));
    end
    
    %% ── 2. DEFINITIONS ───────────────────────────────────────────────────────
    nClasses = size(T.FA_summaryFrac, 2);  % calculates number of columns (= categories)
    classLabels = arrayfun(@(c) sprintf('Cat%d', c), 1:nClasses, 'UniformOutput', false);
    classColors = cell2mat(arrayfun(@(c) getColorClass(c), 1:nClasses, 'UniformOutput', false)');
    
    % legendNames
    if isempty(options.legendNames)
        options.legendNames = classLabels;
    end
    
    %diseaseCats = {'CT','CCM'};
    diseaseCats = cellstr(unique(T.diseaseCat, 'stable'));
    nDis        = numel(diseaseCats);
    
    folders      = unique(T.DataFolder, 'stable');
    nFolders     = numel(folders);
    markerShapes = getPlotMarkers(nFolders);
    if nFolders > numel(markerShapes)
        error('More DataFolders (%d) than available marker shapes (%d).', ...
              nFolders, numel(markerShapes));
    end
    
    % Resolve levels
    if iscategorical(T.extrainfoCat)
        allLevels = cellstr(unique(T.extrainfoCat, 'stable'));
    else
        allLevels = unique(T.extrainfoCat, 'stable');
        if isstring(allLevels), allLevels = cellstr(allLevels); end
    end
    if ~isempty(options.Levels)
        badLevels = options.Levels(~ismember(options.Levels, allLevels));
        if ~isempty(badLevels)
            error('Requested level(s) not found in extrainfoCat: %s', ...
                  strjoin(badLevels, ', '));
        end
        levels = options.Levels;
    else
        levels = allLevels;
    end
    nLevels = numel(levels);
    
    %% ── 3. SELECT METRIC ─────────────────────────────────────────────────────
    if options.Metric == "FRET"
        dataMatrix = T.FA_summaryFRETPx;
        yLabel     = 'Mean FRET per FA';
        titleStr   = 'FRET Distribution by FA Class';
    elseif options.Metric == "Int"
        dataMatrix = T.FA_summaryIntPx;
        yLabel     = 'Mean Intensity per FA (a.u.)';
        titleStr   = 'Intensity Distribution by FA Class';
    else
        dataMatrix = T.FA_summaryFrac;
        yLabel     = 'Fraction of FA';
        titleStr   = 'FA Class Fraction Distribution';
    end

    %% ── 4. PRE-COMPUTE SHARED Y-LIMITS ──────────────────────────────────────
    if options.Metric == "Frac"
        yLims = [0, 1];
    else
        allValsPrecompute = dataMatrix(~isnan(dataMatrix(:)));
        if options.Metric == "Int" && isfinite(options.XLimit)
            allValsPrecompute = allValsPrecompute(allValsPrecompute <= options.XLimit);
        end
        if isempty(allValsPrecompute)
            yLims = [0, 1];
        else
            yPad  = range(allValsPrecompute) * 0.05;
            yLims = [min(allValsPrecompute) - yPad, max(allValsPrecompute) + yPad];
        end
    end
    
    %% ── 5. FIGURE & SUBPLOTS ─────────────────────────────────────────────────
    fig = figure('Color','w','Units','normalized','Position',[0.02 0.05 0.96 0.82]);
    tl = tiledlayout(fig, nDis, nLevels+1, 'TileSpacing','compact','Padding','compact');
    
    axH = gobjects(nDis, nLevels);

    % Data Panels
    for di = 1:nDis
        for ei = 1:nLevels
            %% ── 6. FILTER ROWS FOR THIS PANEL ────────────────────────────
            maskRows = T.diseaseCat == diseaseCats{di} & ...
                       T.extrainfoCat == levels{ei};

            spIdx = (di - 1) * (nLevels+1) + ei;
            ax = nexttile(tl, spIdx);
            axH(di, ei) = ax;
            hold(ax, 'on');

            if sum(maskRows) == 0
                title(ax, sprintf('%s – %s', diseaseCats{di}, levels{ei}), ...
                      'FontSize', 10);
                axis(ax, 'off');
                continue
            end

            %% ── 7. VIOLIN + SCATTER PER CLASS ────────────────────────────
            classData = cell(1, nClasses);

            for c = 1:nClasses
                vals = dataMatrix(maskRows, c);
                vals = vals(~isnan(vals));

                if options.Metric == "Int" && isfinite(options.XLimit)
                    vals = vals(vals <= options.XLimit);
                end

                classData{c} = vals;
                col = classColors(c,:);
                x   = c;

                if numel(vals) >= 2
                    % violin
                    [f, xi] = ksdensity(vals);
                    f       = f / max(f) * 0.35;
                    fill(ax, [x + f, fliplr(x - f)], [xi, fliplr(xi)], col, ...
                        'FaceAlpha', 0.30, ...
                        'EdgeColor', col, ...
                        'EdgeAlpha', 0.5, ...
                        'LineWidth', 0.8);

                    % median line
                    med = median(vals, 'omitnan');
                    plot(ax, [x - 0.18, x + 0.18], [med, med], ...
                        'Color',     col * 0.65, ...
                        'LineWidth', 2.5);

                    % IQR bar
                    iq = quantile(vals, [0.25 0.75]);
                    plot(ax, [x, x], iq, ...
                        'Color',     col * 0.65, ...
                        'LineWidth', 1.5);
                end

                % scatter per DataFolder
                for fi = 1:nFolders
                    folderMask = maskRows & strcmp(T.DataFolder, folders{fi});
                    fVals      = dataMatrix(folderMask, c);
                    fVals      = fVals(~isnan(fVals));

                    if options.Metric == "Int" && isfinite(options.XLimit)
                        fVals = fVals(fVals <= options.XLimit);
                    end
                    if isempty(fVals), continue; end

                    nPts   = numel(fVals);
                    jitter = (rand(nPts,1) - 0.5) * 0.15;
                    scatter(ax, x + jitter, fVals, 22, ...
                        'Marker',            markerShapes{fi}, ...
                        'MarkerFaceColor',   col, ...
                        'MarkerEdgeColor',   col * 0.65, ...
                        'MarkerFaceAlpha',   0.6, ...
                        'MarkerEdgeAlpha',   0.8, ...
                        'LineWidth',         0.6);
                end
            end

            %% ── 8. PANEL FORMATTING ──────────────────────────────────────
            ax.XTick      = 1:nClasses;
            % ax.XTickLabel = classLabels;
            ax.XTickLabel = [];
            ax.XLim       = [0.5, nClasses + 0.5];
            ax.FontSize   = 9;
            ax.LineWidth  = 0.8;
            box(ax, 'off');

            % Set shared y-limits BEFORE drawing brackets so they can expand upward
            ylim(ax, yLims);

            %% ── 9. ONE-WAY ANOVA + SIGNIFICANCE BRACKETS ─────────────────
            if options.ShowStats
                drawAnovaBrackets(ax, classData, options.Alpha);
            end

            if ei == 1
                ylabel(ax, sprintf('%s\n%s', diseaseCats{di}, yLabel), ...
                       'FontSize', 10, 'FontWeight', 'bold');
            end
            if di == 1
                title(ax, levels{ei}, 'FontSize', 11, 'FontWeight', 'bold');
            end

        end  % ei
    end  % di
    
    %% ── 10. FIGURE-LEVEL LABELS & LEGENDS ───────────────────────────────────
    sgtitle(titleStr, 'FontSize', 14, 'FontWeight', 'bold');
    
    % FA Class legend
    axLeg1 = nexttile(tl, nLevels+1);
    cla(axLeg1); axis(axLeg1, 'off');
    hold(axLeg1, 'on');
    xlim(axLeg1, [0 1]); ylim(axLeg1, [0 1]);
    
    text(axLeg1, 0.5, 0.95, 'FA Class', ...
        'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
    
    for c = 1:nClasses
        yPos = 0.80 - (c-1)*0.15;
        patch(axLeg1, ...
            [0.10 0.30 0.30 0.10], ...
            [yPos-0.04, yPos-0.04, yPos+0.04, yPos+0.04], ...
            classColors(c,:), ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.8);
        text(axLeg1, 0.38, yPos, options.legendNames{c}, ...
            'FontSize', 10, 'VerticalAlignment', 'middle');
    end
    
    % Experiment legend
    axLeg2 = nexttile(tl, 2*(nLevels+1));
    cla(axLeg2); axis(axLeg2, 'off');
    hold(axLeg2, 'on');
    xlim(axLeg2, [0 1]); ylim(axLeg2, [0 1]);
    
    text(axLeg2, 0.5, 0.95, 'Experiment', ...
        'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
    
    yStep = min(0.12, 0.85 / nFolders);
    for fi = 1:nFolders
        yPos = 0.85 - (fi-1)*yStep;
        scatter(axLeg2, 0.20, yPos, 36, ...
            'Marker',          markerShapes{fi}, ...
            'MarkerFaceColor', [0.4 0.4 0.4], ...
            'MarkerEdgeColor', [0.2 0.2 0.2]);
        text(axLeg2, 0.32, yPos, char(folders(fi)), ...
            'FontSize', 9, 'VerticalAlignment', 'middle', ...
            'Interpreter', 'none');
    end
end

%% ════════════════════════════════════════════════════════════════════════
%  LOCAL HELPERS
%% ════════════════════════════════════════════════════════════════════════
function drawAnovaBrackets(ax, classData, alpha)
    nClasses = numel(classData);

    allY = [];
    allG = [];
    for c = 1:nClasses
        v = classData{c}(:);
        if ~isempty(v)
            allY = [allY; v];
            allG = [allG; repmat(c, numel(v), 1)];
        end
    end

    validGroups = find(cellfun(@(v) numel(v) >= 2, classData));
    if numel(validGroups) < 2, return; end

    [pAnova, ~, stats] = anova1(allY, allG, 'off');

    if pAnova >= alpha
        yl   = ylim(ax);
        yTop = yl(2);
        text(ax, (nClasses+1)/2, yTop, 'ns (ANOVA)', ...
            'FontSize', 7, 'Color', [0.5 0.5 0.5], ...
            'HorizontalAlignment', 'center');
        return
    end

    [results, ~, ~, ~] = multcompare(stats, ...
        'CType',   'tukey-kramer', ...
        'Alpha',   alpha, ...
        'Display', 'off');

    sigPairs = results( results(:,6) < alpha & ...
                        ismember(results(:,1), validGroups) & ...
                        ismember(results(:,2), validGroups), [1 2 6] );
    if isempty(sigPairs), return; end

    % Sort by span ascending (narrow brackets sit lower)
    spans = abs(sigPairs(:,2) - sigPairs(:,1));
    [~, sortIdx] = sort(spans, 'ascend');
    sigPairs = sigPairs(sortIdx, :);

    % Per-column ceiling: top of actual data (= where the violin tops out)
    yl         = ylim(ax);
    yRange     = yl(2) - yl(1);
    colTop     = zeros(1, nClasses);
    for c = 1:nClasses
        v = classData{c};
        if ~isempty(v)
            colTop(c) = max(v(~isnan(v)), [], 'omitnan');
        else
            colTop(c) = yl(1);
        end
    end
    % Small clearance above the violin
    colTop = colTop + 0.06 * yRange;

    bracketStep = 0.12 * yRange;
    footGap     = 0.03 * yRange;   % gap between a bar and the next bracket's foot

    for k = 1:size(sigPairs, 1)
        g1   = sigPairs(k, 1);
        g2   = sigPairs(k, 2);
        pVal = sigPairs(k, 3);
        cols = g1:g2;

        yBase = max(colTop(cols)) + bracketStep;

        yl = ylim(ax);
        if yBase > yl(2)
            ylim(ax, [yl(1), yBase + bracketStep]);
            yl = ylim(ax);
        end

        yFootA = colTop(g1);
        yFootB = colTop(g2);

        drawSigBracket(ax, g1, g2, yBase, yFootA, yFootB, pVal, alpha);

        % Set ceiling to bar height + footGap so the next foot clears this bar
        colTop(cols) = yBase + footGap;
    end

    % Final small top margin
    yl = ylim(ax);
    ylim(ax, [yl(1), yl(2) + 0.04*(yl(2)-yl(1))]);
end


function drawSigBracket(ax, x1, x2, yBase, yFootA, yFootB, pVal, alpha)
    bracketColor = [0.25 0.25 0.25];
    lineW        = 0.9;

    % Vertical lines from each column's data ceiling up to the horizontal bar
    plot(ax, [x1 x1], [yFootA, yBase], '-', 'Color', bracketColor, 'LineWidth', lineW);
    plot(ax, [x2 x2], [yFootB, yBase], '-', 'Color', bracketColor, 'LineWidth', lineW);
    % Horizontal bar
    plot(ax, [x1 x2], [yBase, yBase],  '-', 'Color', bracketColor, 'LineWidth', lineW);

    yl    = ylim(ax);
    label = pValToStars(pVal, alpha);
    text(ax, (x1+x2)/2, yBase + 0.005*(yl(2)-yl(1)), label, ...
        'FontSize',            8, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'Color',               bracketColor);
end


function stars = pValToStars(p, alpha)
    if     p < alpha / 50,  stars = '***';
    elseif p < alpha / 5,   stars = '**';
    else,                   stars = '*';
    end
end

% function drawAnovaBrackets(ax, classData, alpha)
%     % drawAnovaBrackets  Run one-way ANOVA on classData{1..n}, then draw
%     % significance brackets for all significant pairwise comparisons
%     % (Tukey HSD post-hoc) above the existing violin plots.
%     %
%     % Bracket stacking: pairs are sorted by span width; each bracket is
%     % placed at a y-level that avoids previously drawn brackets.
% 
%     nClasses = numel(classData);
% 
%     % ── Build group vector for anova1 ─────────────────────────────────────
%     allY = [];
%     allG = [];
%     for c = 1:nClasses
%         v = classData{c}(:);
%         if ~isempty(v)
%             allY = [allY; v];           %#ok<AGROW>
%             allG = [allG; repmat(c, numel(v), 1)]; %#ok<AGROW>
%         end
%     end
% 
%     % Need at least 2 non-empty groups with ≥2 observations each
%     validGroups = find(cellfun(@(v) numel(v) >= 2, classData));
%     if numel(validGroups) < 2
%         return
%     end
% 
%     % ── One-way ANOVA (suppress output/figure) ────────────────────────────
%     [pAnova, ~, stats] = anova1(allY, allG, 'off');
% 
%     if pAnova >= alpha
%         % Not significant overall — annotate with 'ns' at top of panel
%         yTop = max(allY) + range(allY)*0.08;
%         text(ax, (nClasses+1)/2, yTop, 'ns (ANOVA)', ...
%             'FontSize', 7, 'Color', [0.5 0.5 0.5], ...
%             'HorizontalAlignment', 'center');
%         return
%     end
% 
%     % ── Tukey HSD post-hoc ────────────────────────────────────────────────
%     [results, ~, ~, ~] = multcompare(stats, ...
%         'CType',   'tukey-kramer', ...
%         'Alpha',   alpha, ...
%         'Display', 'off');
%     % results columns: [group1, group2, lowerCI, diff, upperCI, p-value]
% 
%     % Keep only significant pairs involving valid groups
%     sigPairs = results( results(:,6) < alpha & ...
%                         ismember(results(:,1), validGroups) & ...
%                         ismember(results(:,2), validGroups), [1 2 6] );
% 
%     if isempty(sigPairs)
%         return
%     end
% 
%     % ── Sort pairs by span (narrow first → stack outward) ─────────────────
%     spans = abs(sigPairs(:,2) - sigPairs(:,1));
%     [~, sortIdx] = sort(spans, 'ascend');
%     sigPairs = sigPairs(sortIdx, :);
% 
%     % ── Compute bracket positions ─────────────────────────────────────────
%     yMax = max(allY);
%     yRange = range(allY);
%     bracketStep  = yRange * 0.10;   % vertical gap between bracket levels
%     bracketStart = yMax + yRange * 0.08;
% 
%     % Track occupied y-levels to avoid overlap
%     % Key: pair of x-positions; value: next free y for that x column
%     levelMap = zeros(1, nClasses);  % highest y used per x position
% 
%     for k = 1:size(sigPairs, 1)
%         g1 = sigPairs(k, 1);
%         g2 = sigPairs(k, 2);
%         pVal = sigPairs(k, 3);
% 
%         % y must clear all columns spanned by this bracket
%         cols    = g1:g2;
%         yNeeded = max(levelMap(cols)) + bracketStep;
%         yNeeded = max(yNeeded, bracketStart);
% 
%         drawSigBracket(ax, g1, g2, yNeeded, bracketStep*0.25, pVal, alpha);
% 
%         % Update occupied heights for spanned columns
%         levelMap(cols) = yNeeded + bracketStep * 0.5;
%     end
% end
% 
% 
% function drawSigBracket(ax, x1, x2, yBase, tickH, pVal, alpha)
%     % drawSigBracket  Draw a single significance bracket between x1 and x2
%     % at height yBase, with a label derived from pVal.
% 
%     bracketColor = [0.25 0.25 0.25];
%     lineW        = 0.9;
% 
%     % Bracket lines: left tick | horizontal bar | right tick
%     plot(ax, [x1 x1],   [yBase - tickH, yBase], '-', ...
%         'Color', bracketColor, 'LineWidth', lineW);
%     plot(ax, [x1 x2],   [yBase, yBase],         '-', ...
%         'Color', bracketColor, 'LineWidth', lineW);
%     plot(ax, [x2 x2],   [yBase, yBase - tickH], '-', ...
%         'Color', bracketColor, 'LineWidth', lineW);
% 
%     % Significance label
%     label = pValToStars(pVal, alpha);
% 
%     text(ax, (x1 + x2) / 2, yBase + tickH * 0.3, label, ...
%         'FontSize',           8, ...
%         'HorizontalAlignment','center', ...
%         'VerticalAlignment',  'bottom', ...
%         'Color',              bracketColor);
% end
% 
% 
% function stars = pValToStars(p, alpha)
%     % pValToStars  Convert a p-value to star notation.
%     % Thresholds are relative to alpha (e.g. alpha=0.05):
%     %   p < alpha/50   → ***
%     %   p < alpha/5    → **
%     %   p < alpha      → *
%     %   p >= alpha     → ns
% 
%     if     p < alpha / 50
%         stars = '***';
%     elseif p < alpha / 5
%         stars = '**';
%     else
%         stars = '*';
%     end
% end

% function fig = Plot_FAClassDistribution(T, options)
% % Plot_FAClassDistribution  Violin + scatter of FRET or Intensity per FA
% % class, split by diseaseCat (CT/CCM) × extrainfoCat (rest/R1/R2/1kPa/10kPa).
% % Scatter dot shape varies by DataFolder (experiment).
% % Requires Calc_FAIntClassComparison to have been run first.
% %
% % USAGE:
% %   fig = Plot_FAClassDistribution(T);
% %   fig = Plot_FAClassDistribution(T, Metric='Int');
% %   fig = Plot_FAClassDistribution(T, Metric='Int', XLimit=500);
% %   fig = Plot_FAClassDistribution(T, ShowStats=false);
% % 
% % DEPENDENCIES (must be on MATLAB path)
% %   Filter_IncludeExclude.m
% %   getPlotMarkers.m
% 
% arguments
%     T       table
%     options.Metric  (1,1) string {mustBeMember(options.Metric, ...
%                                  ["FRET","Int", "Frac"])} = "FRET"
%     options.XLimit  (1,1) double = Inf
%     options.Levels      cell          = {}     % ordered subset of extrainfoCat values
%     options.includePatterns cell      = {}
%     options.excludePatterns cell      = {}
%     options.legendNames cell = {}
% end
% 
% % rename
% if isempty(options.legendNames)
%     options.legendNames = {'ThinLine','EllipsoidLine','SmallRound','BigRound'};
% end
% 
% %% Filter rows via include/exclude on SampleFolder
% T = Filter_IncludeExclude(T, options.includePatterns, options.excludePatterns);
% 
% %% ── 1. CHECK PRE-COMPUTED COLUMNS ────────────────────────────────────────
% if options.Metric == "FRET"
%     required = {'FA_summaryFRET'};
% elseif options.Metric == "Int"
%     required = {'FA_summaryInt'};
% else
%     required = {'FA_summaryFrac'};
% end
% missing  = required(~ismember(required, T.Properties.VariableNames));
% if ~isempty(missing)
%     error('Missing columns: %s\nRun Calc_FAIntClassComparison(T) first.', ...
%           strjoin(missing, ', '));
% end
% 
% %% ── 2. DEFINITIONS ───────────────────────────────────────────────────────
% classLabels = {'Cat1','Cat2','Cat3','Cat4'};
% classColors = [
%     getColorClass('Cat1');
%     getColorClass('Cat2');
%     getColorClass('Cat3');
%     getColorClass('Cat4');
% ];
% % classColors = [
% %     0.173, 0.447, 0.698;   % blue   – Cat1
% %     0.302, 0.686, 0.290;   % green  – Cat2
% %     0.890, 0.467, 0.106;   % orange – Cat3
% %     0.506, 0.278, 0.631;   % purple – Cat4
% % ];
% nClasses = numel(classLabels);
% 
% diseaseCats = {'CT','CCM'};
% %extraCats   = {'rest','R1','R2','10kPa','1kPa'};
% nDis        = numel(diseaseCats);   % 2 rows
% %nExtra      = numel(extraCats);     % 5 columns
% 
% % unique experiments → marker shape mapping
% folders      = unique(T.DataFolder, 'stable');
% nFolders     = numel(folders);
% markerShapes = getPlotMarkers(nFolders);
% if nFolders > numel(markerShapes)
%     error('More DataFolders (%d) than available marker shapes (%d).', ...
%           nFolders, numel(markerShapes));
% end
% 
% % RESOLVE LEVELS FROM TABLE
% % Retrieve unique extrainfoCat values in order of appearance
% if iscategorical(T.extrainfoCat)
%     allLevels = cellstr(unique(T.extrainfoCat, 'stable'));
% else
%     allLevels = unique(T.extrainfoCat, 'stable');
%     if isstring(allLevels), allLevels = cellstr(allLevels); end
% end
% % If user supplied Levels, validate and use that order; otherwise use all
% if ~isempty(options.Levels)
%     % validate each requested level exists in table
%     badLevels = options.Levels(~ismember(options.Levels, allLevels));
%     if ~isempty(badLevels)
%         error('Requested level(s) not found in extrainfoCat: %s', ...
%               strjoin(badLevels, ', '));
%     end
%     levels = options.Levels;
% else
%     levels = allLevels;
% end
% nLevels = numel(levels);
% 
% %% ── 3. SELECT METRIC ─────────────────────────────────────────────────────
% if options.Metric == "FRET"
%     dataMatrix = T.FA_summaryFRET;
%     yLabel     = 'Mean FRET per FA';
%     titleStr   = 'FRET Distribution by FA Class';
% elseif options.Metric == "Int"
%     dataMatrix = T.FA_summaryInt;
%     yLabel     = 'Mean Intensity per FA (a.u.)';
%     titleStr   = 'Intensity Distribution by FA Class';
% else
%     dataMatrix = T.FA_summaryFrac;
%     yLabel     = 'Fraction of FA';
%     titleStr   = 'FA Class Fraction Distribution';
% end
% 
% %% ── 4. FIGURE & SUBPLOTS ─────────────────────────────────────────────────
% fig = figure('Color','w','Units','normalized','Position',[0.02 0.05 0.96 0.82]);
% % % shared y-axis limits — compute after collecting all data
% allVals = [];
% tl = tiledlayout(fig, nDis, nLevels+1, 'TileSpacing','compact','Padding','compact');
% 
% % reserve legend tiles (col 6, both rows) - don't create them yet
% axH = gobjects(nDis, nLevels);
% 
% % Data Panels
% for di = 1:nDis
%     for ei = 1:nLevels
% 
%         %% ── 5. FILTER ROWS FOR THIS PANEL ────────────────────────────────
%         maskRows = T.diseaseCat == diseaseCats{di} & ...
%                    T.extrainfoCat == levels{ei};
%         Tsub = T(maskRows, :);
% 
%         spIdx = (di - 1) * (nLevels+1) + ei;
%         ax = nexttile(tl, spIdx);
%         axH(di, ei) = ax;
%         hold(ax, 'on');
% 
%         if height(Tsub) == 0
%             title(ax, sprintf('%s – %s', diseaseCats{di}, levels{ei}), ...
%                   'FontSize', 10);
%             axis(ax, 'off');
%             continue
%         end
% 
%         %% ── 6. VIOLIN + SCATTER PER CLASS ────────────────────────────────
%         for c = 1:nClasses
%             vals = dataMatrix(maskRows, c);
%             vals = vals(~isnan(vals));
% 
%             if options.Metric == "Int" && isfinite(options.XLimit)
%                 vals = vals(vals <= options.XLimit);
%             end
% 
%             col = classColors(c,:);
%             x   = c;
% 
%             allVals = [allVals; vals];
% 
%             if numel(vals) >= 2
%                 % ── violin ────────────────────────────────────────────────
%                 [f, xi] = ksdensity(vals);
%                 f       = f / max(f) * 0.35;
%                 fill(ax, [x + f, fliplr(x - f)], [xi, fliplr(xi)], col, ...
%                     'FaceAlpha', 0.30, ...
%                     'EdgeColor', col, ...
%                     'EdgeAlpha', 0.5, ...
%                     'LineWidth', 0.8);
% 
%                 % ── median line ───────────────────────────────────────────
%                 med = median(vals, 'omitnan');
%                 plot(ax, [x - 0.18, x + 0.18], [med, med], ...
%                     'Color',     col * 0.65, ...
%                     'LineWidth', 2.5);
% 
%                 % ── IQR bar ───────────────────────────────────────────────
%                 iq = quantile(vals, [0.25 0.75]);
%                 plot(ax, [x, x], iq, ...
%                     'Color',     col * 0.65, ...
%                     'LineWidth', 1.5);
%             end
% 
%             % ── scatter per DataFolder (shape-coded) ──────────────────────
%             for fi = 1:nFolders
%                 folderMask = maskRows & strcmp(T.DataFolder, folders{fi});
%                 fVals      = dataMatrix(folderMask, c);
%                 fVals      = fVals(~isnan(fVals));
% 
%                 if options.Metric == "Int" && isfinite(options.XLimit)
%                     fVals = fVals(fVals <= options.XLimit);
%                 end
%                 if isempty(fVals), continue; end
% 
%                 nPts   = numel(fVals);
%                 jitter = (rand(nPts,1) - 0.5) * 0.15;
%                 scatter(ax, x + jitter, fVals, 22, ...
%                     'Marker',            markerShapes{fi}, ...
%                     'MarkerFaceColor',   col, ...
%                     'MarkerEdgeColor',   col * 0.65, ...
%                     'MarkerFaceAlpha',   0.6, ...
%                     'MarkerEdgeAlpha',   0.8, ...
%                     'LineWidth',         0.6);
%             end
%         end
% 
%         %% ── 7. PANEL FORMATTING ──────────────────────────────────────────
%         ax.XTick      = 1:nClasses;
%         ax.XTickLabel = classLabels;
%         ax.XLim       = [0.5, nClasses + 0.5];
%         ax.FontSize   = 9;
%         ax.LineWidth  = 0.8;
%         box(ax, 'off');
% 
%         % row label (disease) on left-most panels only
%         if ei == 1
%             ylabel(ax, sprintf('%s\n%s', diseaseCats{di}, yLabel), ...
%                    'FontSize', 10, 'FontWeight', 'bold');
%         end
% 
%         % column label (extrainfo) on top row only
%         if di == 1
%             title(ax, levels{ei}, 'FontSize', 11, 'FontWeight', 'bold');
%         end
%     end
% end
% 
% %% ── 8. SHARED Y-AXIS LIMITS ──────────────────────────────────────────────
% if ~isempty(allVals)
%     if options.Metric == "Frac"
%         yLims = [0, 1];
%     else
%         yPad  = range(allVals) * 0.05;
%         yLims = [min(allVals) - yPad, max(allVals) + yPad];
%     end
%     for di = 1:nDis
%         for ei = 1:nLevels
%             if isvalid(axH(di,ei)) && strcmp(axH(di,ei).Visible,'on')
%                 ylim(axH(di,ei), yLims);
%             end
%         end
%     end
% end
% 
% %% ── 9. FIGURE-LEVEL LABELS & LEGENDS ────────────────────────────────────
% sgtitle(titleStr, 'FontSize', 14, 'FontWeight', 'bold');
% 
% % ── class colour legend (bottom left) ────────────────────────────────────
% % ── Tile 6: FA Class legend (manual drawing) ─────────────────────────────
% axLeg1 = nexttile(tl, nLevels+1);
% cla(axLeg1); axis(axLeg1, 'off');
% hold(axLeg1, 'on');
% xlim(axLeg1, [0 1]); ylim(axLeg1, [0 1]);
% 
% text(axLeg1, 0.5, 0.95, 'FA Class', ...
%     'FontSize', 10, 'FontWeight', 'bold', ...
%     'HorizontalAlignment', 'center', 'Units', 'normalized');
% 
% for c = 1:nClasses
%     yPos = 0.80 - (c-1)*0.15;
%     patch(axLeg1, ...
%         [0.10 0.30 0.30 0.10], ...
%         [yPos-0.04, yPos-0.04, yPos+0.04, yPos+0.04], ...
%         classColors(c,:), ...
%         'EdgeColor', 'none', ...
%         'FaceAlpha', 0.8);
%     text(axLeg1, 0.38, yPos, options.legendNames{c}, ...
%         'FontSize', 10, 'VerticalAlignment', 'middle');
% end
% 
% % ── Tile 12: Experiment legend (manual drawing) ───────────────────────────
% axLeg2 = nexttile(tl, 2*(nLevels+1));
% cla(axLeg2); axis(axLeg2, 'off');
% hold(axLeg2, 'on');
% xlim(axLeg2, [0 1]); ylim(axLeg2, [0 1]);
% 
% text(axLeg2, 0.5, 0.95, 'Experiment', ...
%     'FontSize', 10, 'FontWeight', 'bold', ...
%     'HorizontalAlignment', 'center', 'Units', 'normalized');
% 
% yStep = min(0.12, 0.85 / nFolders);
% for fi = 1:nFolders
%     yPos = 0.85 - (fi-1)*yStep;
%     scatter(axLeg2, 0.20, yPos, 36, ...
%         'Marker',          markerShapes{fi}, ...
%         'MarkerFaceColor', [0.4 0.4 0.4], ...
%         'MarkerEdgeColor', [0.2 0.2 0.2]);
%     text(axLeg2, 0.32, yPos, char(folders(fi)), ...
%         'FontSize', 9, 'VerticalAlignment', 'middle', ...
%         'Interpreter', 'none');
% end
% end