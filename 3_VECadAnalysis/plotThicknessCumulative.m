function [G, figPanel] = plotThicknessCumulative(T, column, Info, options)
% plotThicknessCumulative
% Plots cumulative distributions (CDF) of Thickness in a 3-row panel
%
% INPUTS (required):
%   T        : table with variables:
%                - T.diseaseCat   (categorical; 'CT','CCM')
%                - T.extrainfoCat (categorical; 'rest','R1','R2','10kPa','1kPa')
%                - T.(column)     (cell of N_i x 1 vectors)
%   column   : string, name of the column to plot
%   xCat     : string, name of the categorical column for levels
%              (used to auto-detect levels when options.levels is empty)
%
% INPUTS (optional, via Name-Value / arguments block):
%   options.binWidth        (1,1) double  = 0.5
%   options.xlimval         (1,1) double  = 10
%   options.includePatterns cell          = {}
%   options.excludePatterns cell          = {}
%   options.levels          cell          = {}   (auto-detected if empty)
%
% OUTPUTS:
%   G        : grouped table with concatenated Thickness vectors
%   figPanel : figure handle of panel figure
%
% EXAMPLE:
%   [G, fig] = plotThicknessCumulative(T, 'Thickness', 'extrainfoCat', ...
%       'binWidth', 0.5, 'xlimval', 10, 'levels', {'rest','R1','R2'});

    arguments
        T                               table
        column                          string
        Info                            struct
        options.xCat                    string = "extrainfoCat"
        options.binWidth        (1,1)   double  = 0.5
        options.xlimval         (1,1)   double  = 10
        options.includePatterns         cell    = {}
        options.excludePatterns         cell    = {}
        options.levels                  cell    = {}
    end

    % Unpack options for convenience
    xCat            = options.xCat;
    binWidth        = options.binWidth;
    xlimval         = options.xlimval;
    includePatterns = options.includePatterns;
    excludePatterns = options.excludePatterns;
    levels          = options.levels;

    % Cast to char so MATLAB built-ins (groupsummary, ismember, etc.) accept them
    column = char(column);
    xCat   = char(xCat);

    % =========================================================================
    % 1. Filter table
    % =========================================================================
    T = Filter_IncludeExclude(T, includePatterns, excludePatterns);

    % =========================================================================
    % 2. Ensure column is a cell of column vectors
    % =========================================================================
    if ~iscell(T.(column))
        T.(column) = num2cell(T.(column), 2);
    end
    T.(column) = cellfun(@(v) v(:), T.(column), 'UniformOutput', false);

    % =========================================================================
    % 3. Build grouped table
    % =========================================================================
    concatFun = @(C) {vertcat(C{:})};
    %G = groupsummary(T, {'diseaseCat', xCat}, concatFun, column);
    G = groupsummary(T, {'diseaseCat', xCat}, concatFun, {column});
    G.Properties.VariableNames(end) = "Thickness";

    % =========================================================================
    % 4. If levels is empty → infer from table using xCat
    % =========================================================================
    tableT = T;  % alias used in the snippet below
    if isempty(levels)
        if ismember(xCat, tableT.Properties.VariableNames)
            levels = reshape(cellstr(string(unique(tableT.(xCat)))), 1, []);
        else
            error("Column '%s' not found in table.", xCat);
        end
    end

    diseaseList = ["CT", "CCM"];

    % =========================================================================
    % 5. Global histogram edges (converted to physical units via pixSize)
    % =========================================================================
    % convert to µm doesn't work yet
    % pixSize = Info.Image_PixSize;
    % allVals   = cellfun(@(v) double(v) * pixSize, G.Thickness, 'UniformOutput', false);
    % allVals   = vertcat(allVals{:});
    % if isempty(allVals)
    %     warning('No Thickness data found.');
    %     return;
    % end
    % globalMin  = floor(min(allVals) / binWidth) * binWidth;
    % globalMax  = ceil(max(allVals)  / binWidth) * binWidth;
    % edges      = globalMin : binWidth : globalMax;
    % if numel(edges) < 2
    %     edges = [globalMin, globalMin + binWidth];
    % end
    % binCenters = edges(1:end-1) + diff(edges) / 2;
    % nBins      = numel(binCenters);
    % 
    % % Convert all stored Thickness vectors to physical units
    % G.Thickness = cellfun(@(v) double(v) * pixSize, G.Thickness, 'UniformOutput', false);

    allVals = vertcat(G.Thickness{:});
    if isempty(allVals)
        warning('No Thickness data found.');
        return;
    end
    globalMin = floor(min(allVals) / binWidth) * binWidth;
    globalMax = ceil(max(allVals)  / binWidth) * binWidth;
    edges     = globalMin : binWidth : globalMax;
    if numel(edges) < 2
        edges = [globalMin, globalMin + binWidth];
    end
    binCenters = edges(1:end-1) + diff(edges) / 2;
    nBins      = numel(binCenters);

    % =========================================================================
    % 6. Helper: compute CDF + SD band per group
    %    Returns cumProb (1 x nBins) and sdBand (1 x nBins).
    %    SD is computed by bootstrap-resampling subjects in a group.
    % =========================================================================
    function [cumProb, sdBand] = computeCDF(vals, edgesIn)
        if isempty(vals)
            cumProb = nan(1, numel(edgesIn)-1);
            sdBand  = nan(1, numel(edgesIn)-1);
            return;
        end
        counts  = histcounts(double(vals), edgesIn, 'Normalization', 'probability');
        cumProb = cumsum(counts) * 100;

        % SD via simple ±1 SD across bins using a Bernoulli approximation:
        % At each bin centre x, p = cumProb/100.  Var(p_hat) ≈ p(1-p)/N.
        N       = numel(vals);
        p       = cumProb / 100;
        sdBand  = sqrt(p .* (1 - p) / N) * 100;   % convert back to %
    end

    % =========================================================================
    % 7. Panel figure  (3 rows × 5 columns)
    % =========================================================================
    nLevels  = numel(levels);
    nCols    = max(nLevels, 5);          % keep at least 5 columns
    figPanel = figure('Color', 'w');
    tLayout  = tiledlayout(3, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    % ------------------------------------------------------------------
    % ROW 1: CT vs CCM overlay per extra-level
    % ------------------------------------------------------------------
    for e = 1 : nLevels
        nexttile; hold on; grid on;
        title(levels{e});

        cdfStore   = nan(numel(diseaseList), nBins);
        sdStore    = nan(numel(diseaseList), nBins);
        validLines = false(1, numel(diseaseList));

        for d = 1 : numel(diseaseList)
            mask = string(G.diseaseCat) == diseaseList(d) & ...
                   string(G.(xCat))     == levels{e};
            if any(mask)
                vals = G.Thickness{mask};
                if isempty(vals), continue; end

                [cumProb, sdBand] = computeCDF(vals, edges);
                cdfStore(d, :)    = cumProb;
                sdStore(d, :)     = sdBand;
                validLines(d)     = true;

                clr = 'b';
                if diseaseList(d) == "CCM", clr = 'r'; end

                % SD shaded band
                fill([binCenters, fliplr(binCenters)], ...
                     [cumProb + sdBand, fliplr(cumProb - sdBand)], ...
                     clr, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');

                plot(binCenters, cumProb, '-', 'Color', clr, ...
                     'LineWidth', 1.8, 'DisplayName', char(diseaseList(d)));
            end
        end

        % Max difference between CT and CCM
        % if all(validLines)
        %     diffVec = abs(cdfStore(1,:) - cdfStore(2,:));
        %     [maxDiff, maxIdx] = max(diffVec);
        %     xMaxDiff = binCenters(maxIdx);
        %     yMaxDiff = mean(cdfStore(:, maxIdx));
        %     plot([xMaxDiff xMaxDiff], cdfStore(:, maxIdx)', '--k', ...
        %          'LineWidth', 1, 'HandleVisibility', 'off');
        %     text(xMaxDiff, yMaxDiff, sprintf('  \\DeltaMax=%.1f%%', maxDiff), ...
        %          'FontSize', 7, 'Color', 'k', 'VerticalAlignment', 'middle');
        % end
        % Max horizontal difference between CT and CCM CDFs:
        % For each % level p, find x_CT(p) and x_CCM(p) by inverting each CDF,
        % then report max |x_CT - x_CCM| in thickness units.
        if all(validLines)
            pLevels = 1 : 99;
            % Invert CDF: interpolate x given p — requires unique p values
            cdf1 = cdfStore(1,:);  xc1 = binCenters;
            cdf2 = cdfStore(2,:);  xc2 = binCenters;
            [cdf1u, ia1] = unique(cdf1, 'last');  xc1u = xc1(ia1);
            [cdf2u, ia2] = unique(cdf2, 'last');  xc2u = xc2(ia2);
            xCT  = interp1(cdf1u, xc1u, pLevels, 'linear', nan);
            xCCM = interp1(cdf2u, xc2u, pLevels, 'linear', nan);
            hDiff = abs(xCT - xCCM);
            [maxHDiff, maxPIdx] = max(hDiff);
            pStar = pLevels(maxPIdx);
            xA    = xCT(maxPIdx);
            xB    = xCCM(maxPIdx);
            plot([xA xB], [pStar pStar], '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
            plot([xA xA], [pStar pStar], 'ok', 'MarkerSize', 3, 'HandleVisibility', 'off');
            plot([xB xB], [pStar pStar], 'ok', 'MarkerSize', 3, 'HandleVisibility', 'off');
            text(mean([xA xB]), pStar, sprintf('  \\Deltax=%.2f', maxHDiff), ...
                 'FontSize', 7, 'Color', 'k', 'VerticalAlignment', 'bottom');
        end

        xlabel('Thickness (px)'); ylabel('Cumulative %');
        ylim([0 100]); xlim([0 xlimval]);
        legend({'CT', 'CCM'}, 'Location', 'best');
    end

    % Fill any unused tiles in row 1
    for k = (nLevels+1) : nCols
        nexttile; axis off;
    end

    % ------------------------------------------------------------------
    % Helper: plot a combined-levels row (rows 2 & 3)
    % ------------------------------------------------------------------
    function plotCombinedRow(disease, groupList)
        % groupList : cell array of cell arrays of level strings
        for g = 1 : numel(groupList)
            nexttile; hold on; grid on;
            groupLevels = groupList{g};

            % cdfRest  = [];
            % sdRest   = [];
            % restDone = false;
            cdfRest      = [];
            sdRest       = [];
            restDone     = false;
            allCDFs      = {};   % non-rest CDFs
            nonRestLabels = {};  % matching level names

            allCDFs  = {};  % store non-rest CDFs for max-diff computation

            for lev = 1 : numel(groupLevels)
                mask = string(G.diseaseCat) == disease & ...
                       string(G.(xCat))     == groupLevels{lev};
                if ~any(mask), continue; end
                vals = G.Thickness{mask};
                if isempty(vals), continue; end

                [cumProb, sdBand] = computeCDF(vals, edges);
                clr = getColor(groupLevels{lev});

                % SD shaded band
                fill([binCenters, fliplr(binCenters)], ...
                     [cumProb + sdBand, fliplr(cumProb - sdBand)], ...
                     clr, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');

                plot(binCenters, cumProb, '-', 'Color', clr, ...
                     'LineWidth', 1.8, 'DisplayName', groupLevels{lev});

                if strcmp(groupLevels{lev}, 'rest')
                    cdfRest  = cumProb;
                    sdRest   = sdBand;
                    restDone = true;
                else
                    allCDFs{end+1}       = cumProb;
                    nonRestLabels{end+1} = groupLevels{lev};
                end
            end

            % Max difference: rest vs each other line
            % if restDone && ~isempty(allCDFs)
            %     maxDiffAll = 0;
            %     for ci = 1 : numel(allCDFs)
            %         diffVec = abs(cdfRest - allCDFs{ci});
            %         [md, mi] = max(diffVec);
            %         if md > maxDiffAll
            %             maxDiffAll = md;
            %             maxIdxAll  = mi;
            %             otherCDF   = allCDFs{ci};
            %         end
            %     end
            %     xMaxDiff = binCenters(maxIdxAll);
            %     yVals    = [cdfRest(maxIdxAll), otherCDF(maxIdxAll)];
            %     plot([xMaxDiff xMaxDiff], yVals, '--k', ...
            %          'LineWidth', 1, 'HandleVisibility', 'off');
            %     text(xMaxDiff, mean(yVals), ...
            %          sprintf('  \\DeltaMax=%.1f%%', maxDiffAll), ...
            %          'FontSize', 7, 'Color', 'k', 'VerticalAlignment', 'middle');
            % end
            % Max horizontal diff: rest vs each other line individually
            if restDone && ~isempty(allCDFs)
                pLevels = 1 : 99;
                [cdfRu, iru] = unique(cdfRest, 'last');
                xRestU = binCenters(iru);
                xRest  = interp1(cdfRu, xRestU, pLevels, 'linear', nan);
                vOffsets = linspace(8, -8, numel(allCDFs));
                for ci = 1 : numel(allCDFs)
                    [cdfOu, iou] = unique(allCDFs{ci}, 'last');
                    xOtherU = binCenters(iou);
                    xOther  = interp1(cdfOu, xOtherU, pLevels, 'linear', nan);
                    hDiff   = abs(xRest - xOther);
                    [maxHDiff, maxPIdx] = max(hDiff);
                    pStar   = pLevels(maxPIdx);
                    xA      = xRest(maxPIdx);
                    xB      = xOther(maxPIdx);
                    clrPair = getColor(nonRestLabels{ci});
                    plot([xA xB], [pStar pStar], '--', 'Color', clrPair, ...
                         'LineWidth', 1, 'HandleVisibility', 'off');
                    plot([xA xA; xB xB], [pStar pStar; pStar pStar], 'o', ...
                         'Color', clrPair, 'MarkerSize', 3, 'HandleVisibility', 'off');
                    text(mean([xA xB]), pStar + vOffsets(ci), ...
                         sprintf('  \\Deltax=%.2f\n  (rest vs %s)', maxHDiff, nonRestLabels{ci}), ...
                         'FontSize', 7, 'Color', clrPair, 'VerticalAlignment', 'bottom');
                end
            end

            xlabel('Thickness'); ylabel('Cumulative %');
            ylim([0 100]); xlim([0 xlimval]);
            title(sprintf('%s: %s', disease, strjoin(groupLevels, '+')));
            legend('Location', 'best');
        end

        % Fill unused tiles in this row
        usedTiles = numel(groupList);
        for k = (usedTiles+1) : nCols
            nexttile; axis off;
        end
    end

    % ------------------------------------------------------------------
    % ROW 2: CT combined across multi-extra-levels
    % ------------------------------------------------------------------
    row2Groups = { {'rest','R1','R2'}, {'rest','10kPa','1kPa'} };
    plotCombinedRow("CT", row2Groups);

    % ------------------------------------------------------------------
    % ROW 3: CCM combined across multi-extra-levels
    % ------------------------------------------------------------------
    row3Groups = { {'rest','R1','R2'}, {'rest','10kPa','1kPa'} };
    plotCombinedRow("CCM", row3Groups);

    sgtitle('Cumulative Thickness distributions (%)');

end

% =========================================================================
%  Local helper: getColor
% =========================================================================
function c = getColor(extraLevel)
    switch char(extraLevel)
        case 'rest',   c = [0.0 0.0 0.0];      % black
        case 'R1',     c = [0.6 0.7 1.0];      % light purple
        case 'R2',     c = [0.5 0.1 0.7];      % dark purple
        case '10kPa',  c = [0.6 1.0 0.6];      % light green
        case '1kPa',   c = [0.0 0.5 0.0];      % dark green
        otherwise,     c = [0.5 0.5 0.5];      % grey
    end
end

%% Old Code
% function [G, figPanel] = plotThicknessCumulative(T, column, binWidth, xlimval, includePatterns, excludePatterns)
% % plotThicknessCumulative
% % Plots cumulative distributions (CDF) of Thickness in a 3-row panel
% %
% % INPUTS:
% %   T        : table with variables:
% %                - T.diseaseCat   (categorical; 'CT','CCM')
% %                - T.extrainfoCat (categorical; 'rest','R1','R2','10kPa','1kPa')
% %                - T.Thickness   (cell of N_i x 1 vectors)
% %   binWidth : histogram bin width
% %   xlimval  : x-axis limit
% %
% % OUTPUTS:
% %   G        : grouped table with concatenated Thickness vectors
% %   figPanel : figure handle of panel figure
% 
%     T = Filter_IncludeExclude(T, includePatterns, excludePatterns);
% 
%     if nargin < 2
%         error('Usage: [G, figPanel] = plotThicknessCumulative(T, binWidth, xlimval)');
%     end
% 
%     % Ensure Thickness is cell of column vectors
%     if ~iscell(T.(column))
%         T.(column) = num2cell(T.(column), 2);
%     end
%     T.(column) = cellfun(@(v) v(:), T.(column), 'UniformOutput', false);
% 
%     % Build grouped table
%     concatFun = @(C) {vertcat(C{:})};
%     G = groupsummary(T, {'diseaseCat','extrainfoCat'}, concatFun, column);
%     G.Properties.VariableNames(end) = "Thickness";
% 
%     diseaseList = ["CT","CCM"];
%     infoList    = ["rest","R1","R2","10kPa","1kPa"];
% 
%     % Global histogram edges
%     allVals = vertcat(G.Thickness{:});
%     if isempty(allVals)
%         warning('No Thickness data found.');
%         return;
%     end
%     globalMin = floor(min(allVals)/binWidth)*binWidth;
%     globalMax = ceil(max(allVals)/binWidth)*binWidth;
%     edges = globalMin:binWidth:globalMax;
%     if numel(edges)<2
%         edges = [globalMin, globalMin+binWidth];
%     end
%     binCenters = edges(1:end-1) + diff(edges)/2;
% 
%     %% Panel figure 3x5
%     figPanel = figure('Color','w');
%     tLayout = tiledlayout(3,5,'TileSpacing','compact','Padding','compact');
% 
%     % -------- Row1: CT vs CCM overlay per extraLevel --------
%     for e = 1:numel(infoList)
%         nexttile; hold on; grid on;
%         title(infoList(e));
%         for d = 1:numel(diseaseList)
%             mask = string(G.diseaseCat)==diseaseList(d) & string(G.extrainfoCat)==infoList(e);
%             if any(mask)
%                 vals = G.Thickness{mask};
%                 if isempty(vals), continue; end
%                 counts = histcounts(double(vals), edges, 'Normalization','probability');
%                 cumProb = cumsum(counts) * 100; % percentage
%                 color = 'b'; if diseaseList(d)=="CCM", color='r'; end
%                 plot(binCenters, cumProb, '-', 'Color', color, 'LineWidth',1.8,'DisplayName',string(diseaseList(d)));
%             end
%         end
%         xlabel('Thickness'); ylabel('Cumulative %'); ylim([0 100]); xlim([0 xlimval]);
%         legend({'CT','CCM'}, 'Location','best');
%     end
% 
%     % -------- Row2: CT combined across multi-extraLevels --------
%     row2Groups = { ["rest","R1","R2"], ["rest","10kPa","1kPa"] };
%     for g = 1:numel(row2Groups)
%         nexttile; hold on; grid on;
%         groupLevels = row2Groups{g};
%         for lev = 1:numel(groupLevels)
%             mask = string(G.diseaseCat)=="CT" & string(G.extrainfoCat)==groupLevels(lev);
%             if any(mask)
%                 vals = G.Thickness{mask};
%                 if isempty(vals), continue; end
%                 counts = histcounts(double(vals), edges, 'Normalization','probability');
%                 cumProb = cumsum(counts) * 100;
%                 plot(binCenters, cumProb, '-', 'Color', getColor(groupLevels(lev)), 'LineWidth',1.8, 'DisplayName', groupLevels(lev));
%             end
%         end
%         xlabel('Thickness'); ylabel('Cumulative %'); ylim([0 100]); xlim([0 xlimval]);
%         title(sprintf('CT combined: %s', strjoin(groupLevels,'+')));
%         legend('Location','best');
%     end
%     % Fill empty tiles in row2 (3,4,5) to prevent CCM shifting into row2
%     for k = 3:5
%         nexttile; axis off;
%     end
% 
%     % -------- Row3: CCM combined across multi-extraLevels --------
%     row3Groups = { ["rest","R1","R2"], ["rest","10kPa","1kPa"] };
%     for g = 1:numel(row3Groups)
%         nexttile; hold on; grid on;
%         groupLevels = row3Groups{g};
%         for lev = 1:numel(groupLevels)
%             mask = string(G.diseaseCat)=="CCM" & string(G.extrainfoCat)==groupLevels(lev);
%             if any(mask)
%                 vals = G.Thickness{mask};
%                 if isempty(vals), continue; end
%                 counts = histcounts(double(vals), edges, 'Normalization','probability');
%                 cumProb = cumsum(counts) * 100;
%                 plot(binCenters, cumProb, '-', 'Color', getColor(groupLevels(lev)), 'LineWidth',1.8, 'DisplayName', groupLevels(lev));
%             end
%         end
%         xlabel('Thickness'); ylabel('Cumulative %'); ylim([0 100]); xlim([0 xlimval]);
%         title(sprintf('CCM combined: %s', strjoin(groupLevels,'+')));
%         legend('Location','best');
%     end
% 
%     sgtitle('Cumulative Thickness distributions (%)');
% 
% end
% 
% %% Helper: getColor function
% function c = getColor(extraLevel)
%     switch char(extraLevel)  % ensure char
%         case 'rest',  c = [0 0 0];       % black
%         case 'R1',    c = [0.6 0.7 1];   % light purple
%         case 'R2',    c = [0.5 0.1 0.7]; % dark purple
%         case '10kPa',  c = [0.6 1 0.6];  % light green
%         case '1kPa',   c = [0 0.5 0];    % dark green
%         otherwise,   c = [0.5 0.5 0.5]; % gray for safety
%     end
% end
