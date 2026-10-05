% keep = T;
% 
% 
% % Ensure Thickness is a cell array of column vectors
% if ~iscell(T.Thickness)
%     T.Thickness = num2cell(T.Thickness, 2);
% end
% T.Thickness = cellfun(@(v) v(:), T.Thickness, 'UniformOutput', false);
% 
% % Custom concatenator returns a single cell (so groupsummary can store it)
% concatFun = @(C) {vertcat(C{:})};
% 
% % One row per (diseaseCat, extrainfoCat) with concatenated vector
% G = groupsummary(T, {'diseaseCat','extrainfoCat'}, concatFun, 'Thickness');
% 
% % Optional: make the output variable name nicer
% G.Properties.VariableNames(end) = "Thickness";

%%
% function plotThickness(G, extraLevel)
% % G: table from groupsummary with variables: diseaseCat, extrainfoCat, Thickness (cell of vectors)
% % extraLevel: categorical scalar matching a level in G.extrainfoCat (e.g., categorical("rest"))
% 
%     % Select the two series for this extrainfo level
%     isExtra = G.extrainfoCat == extraLevel;
% 
%     % Fetch CT and CCM rows (may be missing in some datasets)
%     iCT  = find(isExtra & G.diseaseCat == 'CT', 1);
%     iCCM = find(isExtra & G.diseaseCat == 'CCM', 1);
% 
%     % Prepare figure
%     figure('Color','w'); hold on; grid on;
%     set(gca, 'FontSize', 12);
% 
%     % Plot CT if present
%     h = gobjects(1,2); % for legend handles
%     if ~isempty(iCT)
%         yCT = G.Thickness{iCT};
%         xCT = (1:numel(yCT))';
%         h(1) = plot(xCT, double(yCT), 'b-', 'LineWidth', 1.5, 'DisplayName','CT');
%     else
%         warning('CT series missing for extrainfoCat = %s', string(extraLevel));
%     end
% 
%     % Plot CCM if present
%     if ~isempty(iCCM)
%         yCCM = G.Thickness{iCCM};
%         xCCM = (1:numel(yCCM))';
%         h(2) = plot(xCCM, double(yCCM), 'r-', 'LineWidth', 1.5, 'DisplayName','CCM');
%     else
%         warning('CCM series missing for extrainfoCat = %s', string(extraLevel));
%     end
% 
%     % Labels, legend, title
%     xlabel('Index', 'FontSize', 12);
%     ylabel('Prevalence', 'FontSize', 12);
%     title(sprintf('Prevalence vs Index — extrainfo: %s', string(extraLevel)));
% 
%     % Legend only for the plotted lines
%     plotted = isgraphics(h);
%     if any(plotted)
%         legend(h(plotted), 'Location', 'best');
%     end
% end
% 
% 

%%


%% WORKS!
% function [G, fig] = plotThickness(T, extraLevel, binWidth, xlimval)
% % plotThickness
% %  - Builds grouped table G combining Nx1 vectors per (diseaseCat, extrainfoCat)
% %  - Plots CT (blue) vs CCM (red) for one extrainfoCat with histogram binning (bin width = binWidth)
% %
% % INPUTS
% %   T          : table with variables:
% %                   - T.diseaseCat   (categorical; includes 'CT','CCM')
% %                   - T.extrainfoCat (categorical; e.g., 'rest','R1','R2','10kPa','1kP', ...)
% %                   - T.Thickness    (cell of N_i x 1 vectors; often single)
% %   extraLevel : scalar text (char/string) or scalar categorical indicating which extrainfoCat to plot
% %   binWidth   : positive numeric scalar (e.g., 2)
% %
% % OUTPUTS
% %   G   : grouped table with one row per (diseaseCat, extrainfoCat) and concatenated Thickness vector
% %   fig : figure handle
% 
%     % --------- Basic input checks ----------
%     if nargin < 3
%         error('Usage: [G, fig] = plotThickness(T, extraLevel, binWidth)');
%     end
%     validateattributes(T, {'table'}, {'nonempty'}, mfilename, 'T', 1);
%     validateattributes(binWidth, {'numeric'}, {'scalar','positive','finite'}, mfilename, 'binWidth', 3);
% 
%     % Normalize extraLevel to scalar text (avoid ordinal/non-ordinal pitfalls)
%     extraTxt = normalizeExtraAsText(extraLevel);
% 
%     % --------- Ensure Thickness column is a cell array of column vectors ----------
%     if ~iscell(T.Thickness)
%         % If rows contain vectors in a numeric matrix, convert row-wise into cells
%         T.Thickness = num2cell(T.Thickness, 2);
%     end
%     T.Thickness = cellfun(@(v) v(:), T.Thickness, 'UniformOutput', false);
%     % If you want to enforce single dtype across the board, uncomment:
%     % T.Thickness = cellfun(@single, T.Thickness, 'UniformOutput', false);
% 
%     % --------- Build grouped table G: one row per (diseaseCat, extrainfoCat) ----------
%     concatFun = @(C) {vertcat(C{:})}; % groupsummary expects function to return a scalar value; wrap in cell
%     G = groupsummary(T, {'diseaseCat','extrainfoCat'}, concatFun, 'Thickness');
%     G.Properties.VariableNames(end) = "Thickness"; % rename the summary variable
% 
%     % --------- Extract CT / CCM for the requested extraLevel ----------
%     % Compare using strings to avoid ordinal vs non-ordinal categorical issues
%     isExtra = string(G.extrainfoCat) == string(extraTxt);
% 
%     % We will also select disease using strings to avoid any categorical quirks
%     isCT  = string(G.diseaseCat) == "CT";
%     isCCM = string(G.diseaseCat) == "CCM";
% 
%     iCT  = find(isExtra & isCT,  1);
%     iCCM = find(isExtra & isCCM, 1);
% 
%     if isempty(iCT) && isempty(iCCM)
%         % Show available levels to help debugging
%         warning('No CT or CCM data found for extrainfoCat = "%s". Available: %s', ...
%             extraTxt, strjoin(unique(string(G.extrainfoCat))', ', '));
%         error('No matching data for requested extrainfoCat.');
%     end
% 
%     yCT  = []; if ~isempty(iCT),  yCT  = G.Thickness{iCT};  end
%     yCCM = []; if ~isempty(iCCM), yCCM = G.Thickness{iCCM}; end
% 
%     % --------- Shared histogram edges using binWidth ----------
%     %binWidth = 0.5; % HERE BIN!
%     mins = []; maxs = [];
%     if ~isempty(yCT),  mins(end+1) = min(yCT);  maxs(end+1) = max(yCT);  end
%     if ~isempty(yCCM), mins(end+1) = min(yCCM); maxs(end+1) = max(yCCM); end
% 
%     globalMin = floor(min(mins)/binWidth)*binWidth;
%     globalMax = ceil(max(maxs)/binWidth)*binWidth;
%     if ~isfinite(globalMin) || ~isfinite(globalMax) || globalMin == globalMax
%         % Degenerate case when empty or identical values; expand by one bin on each side
%         if ~isfinite(globalMin), globalMin = 0; end
%         if ~isfinite(globalMax), globalMax = 0; end
%         globalMin = globalMin - binWidth;
%         globalMax = globalMax + binWidth;
%     end
%     edges = globalMin:binWidth:globalMax;
% 
%     % Ensure at least two edges
%     if numel(edges) < 2
%         edges = [globalMin, globalMin+binWidth];
%     end
% 
%     % --------- Compute prevalence (probability per bin) ----------
%     ctProb  = [];
%     ccmProb = [];
%     if ~isempty(yCT)
%         ctProb = histcounts(double(yCT), edges, 'Normalization','probability');
%     end
%     if ~isempty(yCCM)
%         ccmProb = histcounts(double(yCCM), edges, 'Normalization','probability');
%     end
%     binCenters = edges(1:end-1) + diff(edges)/2;
% 
%     % --------- Plot ----------
%     fig = figure('Color','w'); hold on; grid on;
%     set(gca, 'FontSize', 12);
% 
%     h = gobjects(1,2);
%     if ~isempty(ctProb)
%         h(1) = plot(binCenters, ctProb, 'b-', 'LineWidth', 1.8, 'DisplayName','CT');
%         % Alternative (stairs): stairs(edges, [ctProb, 0], 'b-', 'LineWidth', 1.8, 'DisplayName','CT');
%     end
%     if ~isempty(ccmProb)
%         h(2) = plot(binCenters, ccmProb, 'r-', 'LineWidth', 1.8, 'DisplayName','CCM');
%         % Alternative (stairs): stairs(edges, [ccmProb, 0], 'r-', 'LineWidth', 1.8, 'DisplayName','CCM');
%     end
% 
%     set(gca, 'YScale', 'log'); % set y-axis to log (comment out if non needed)
% 
%     xlabel(sprintf('Value (bin width = %g)', binWidth));
%     ylabel('Prevalence (probability)');
%     title(sprintf('CT vs CCM — extrainfo: %s', extraTxt));
% 
%     xlim([0, xlimval]);
% 
%     plotted = isgraphics(h);
%     if any(plotted)
%         legend(h(plotted), 'Location','best');
%     end
% 
%     % % Give a bit of headroom on Y
%     % if ~isempty([ctProb, ccmProb])
%     %     ymax = max([ctProb(:); ccmProb(:)]);
%     %     if isfinite(ymax) && ymax > 0
%     %         ylim([0, ymax*1.1]);
%     %     end
%     % end
% end
% 
% % ===== Helper: normalize extraLevel to scalar text =====
% function txt = normalizeExtraAsText(extraLevelIn)
%     if iscategorical(extraLevelIn)
%         if ~isscalar(extraLevelIn)
%             error('extraLevel must be a scalar categorical.');
%         end
%         txt = string(extraLevelIn);
%     elseif ischar(extraLevelIn) || (isstring(extraLevelIn) && isscalar(extraLevelIn))
%         txt = string(extraLevelIn);
%     else
%         error('extraLevel must be a single char, string, or scalar categorical (e.g., ''rest'').');
%     end
% end


%% Also works - update
% function [G, figPanel, figOverlay] = plotThickness(T, binWidth, xlimval)
% % plotThicknessOverlay
% % Plots Thickness histograms:
% % 1) Panel: each extraLevel (1 row, 5 columns) with CT (blue) vs CCM (red)
% % 2) Overlay: all CTs and all CCMs combined in one figure
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
% %   G           : grouped table with concatenated Thickness vectors
% %   figPanel    : figure handle of panel figure
% %   figOverlay  : figure handle of overlay figure
% 
%     if nargin < 2
%         error('Usage: [G, figPanel, figOverlay] = plotThicknessOverlay(T, binWidth, xlimval)');
%     end
% 
%     % Ensure Thickness is a cell of column vectors
%     if ~iscell(T.Thickness)
%         T.Thickness = num2cell(T.Thickness, 2);
%     end
%     T.Thickness = cellfun(@(v) v(:), T.Thickness, 'UniformOutput', false);
% 
%     % -------- Build grouped table ----------
%     concatFun = @(C) {vertcat(C{:})};
%     G = groupsummary(T, {'diseaseCat','extrainfoCat'}, concatFun, 'Thickness');
%     G.Properties.VariableNames(end) = "Thickness";
% 
%     diseaseList = ["CT","CCM"];
%     infoList    = ["rest","R1","R2","10kPa","1kPa"];
% 
%     % -------- Global edges for histograms ----------
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
%     epsVal = 1e-6; % for log scale zeros
% 
%     %% ----------- Panel figure: 1 row, overlay CT vs CCM ----------
%     figPanel = figure('Color','w');
%     tiledlayout(1,5,'TileSpacing','compact','Padding','compact');
% 
%     for e = 1:numel(infoList)
%         nexttile;
%         hold on; grid on;
%         title(infoList(e));
% 
%         for d = 1:numel(diseaseList)
%             mask = string(G.diseaseCat)==diseaseList(d) & string(G.extrainfoCat)==infoList(e);
%             if any(mask)
%                 vals = G.Thickness{mask};
%                 if isempty(vals)
%                     continue;
%                 end
%                 prob = histcounts(double(vals), edges, 'Normalization','probability');
%                 prob(prob==0) = epsVal; % log scale
%                 color = 'b'; if diseaseList(d)=="CCM", color='r'; end
%                 plot(binCenters, prob, '-', 'Color', color, 'LineWidth', 1.8, 'DisplayName', string(diseaseList(d)));
%             end
%         end
% 
%         set(gca,'YScale','log','XLim',[0 xlimval]);
%         xlabel('Thickness'); ylabel('Prevalence');
%         legend({'CT','CCM'}, 'Location','best');
%     end
%     sgtitle('Thickness distributions per extraLevel (CT vs CCM)');
% 
%     %% ----------- Overlay figure: all CTs combined, all CCMs combined ----------
%     figOverlay = figure('Color','w'); hold on; grid on;
%     colors = {'b','r'};
%     labels = {'CT','CCM'};
% 
%     for d = 1:numel(diseaseList)
%         mask = string(G.diseaseCat)==diseaseList(d);
%         valsAll = vertcat(G.Thickness{mask});
%         if ~isempty(valsAll)
%             prob = histcounts(double(valsAll), edges, 'Normalization','probability');
%             prob(prob==0) = epsVal;
%             plot(binCenters, prob, '-', 'Color', colors{d}, 'LineWidth', 2, 'DisplayName', labels{d});
%         end
%     end
%     set(gca,'YScale','log','XLim',[0 xlimval]);
%     xlabel('Thickness'); ylabel('Prevalence');
%     legend('Location','best');
%     title('Overlay of all CTs vs all CCMs');
% end
%%
function [G, figPanel] = plotThickness(T, column, binWidth, xlimval, includePatterns, excludePatterns)
% plotThicknessMultiPanel
% 3-row panel:
% Row1: CT+CCM overlay per extraLevel
% Row2: CT combined across multi-extraLevels
% Row3: CCM combined across multi-extraLevels
% Explicit empty panels used to force layout

    T = Filter_IncludeExclude(T, includePatterns, excludePatterns);

    if nargin < 2
        error('Usage: [G, figPanel] = plotThicknessMultiPanel(T, binWidth, xlimval)');
    end

    % Ensure Thickness is cell of column vectors
    if ~iscell(T.(column))
        T.(column) = num2cell(T.(column), 2);
    end
    T.(column) = cellfun(@(v) v(:), T.(column), 'UniformOutput', false);

    % Build grouped table
    concatFun = @(C) {vertcat(C{:})};
    G = groupsummary(T, {'diseaseCat','extrainfoCat'}, concatFun, column);
    G.Properties.VariableNames(end) = "Thickness";

    diseaseList = ["CT","CCM"];
    infoList    = ["rest","R1","R2","10kPa","1kPa"];

    % Global histogram edges
    allVals = vertcat(G.Thickness{:});
    if isempty(allVals)
        warning('No Thickness data found.');
        return;
    end
    globalMin = floor(min(allVals)/binWidth)*binWidth;
    globalMax = ceil(max(allVals)/binWidth)*binWidth;
    edges = globalMin:binWidth:globalMax;
    if numel(edges)<2
        edges = [globalMin, globalMin+binWidth];
    end
    binCenters = edges(1:end-1) + diff(edges)/2;
    epsVal = 1e-6;

    % Color map for extraLevels
    % colorMap = containers.Map( ...
    %     {"rest","R1","R2","10kPa","1kPa"}, ...
    %     {[0 0 0], [0.6 0.8 1], [0 0 0.8], [1 1 0], [1 0.5 0]});

    %% Panel figure 3x5
    figPanel = figure('Color','w');
    tLayout = tiledlayout(3,5,'TileSpacing','compact','Padding','compact');

    % -------- Row1: CT+CCM per extraLevel --------
    for e = 1:numel(infoList)
        nexttile;
        hold on; grid on; title(infoList(e));

        for d = 1:numel(diseaseList)
            mask = string(G.diseaseCat)==diseaseList(d) & string(G.extrainfoCat)==infoList(e);
            if any(mask)
                vals = G.Thickness{mask};
                if isempty(vals), continue; end
                prob = histcounts(double(vals), edges, 'Normalization','probability');
                prob(prob==0) = epsVal;
                color = 'b'; if diseaseList(d)=="CCM", color='r'; end
                plot(binCenters, prob, '-', 'Color', color, 'LineWidth', 1.8, 'DisplayName', string(diseaseList(d)));
            end
        end
        set(gca,'YScale','log','XLim',[0 xlimval]);
        xlabel('Thickness'); ylabel('Prevalence');
        legend({'CT','CCM'}, 'Location','best');
    end

    % -------- Row2: CT multi-extraLevel --------
    row2Groups = { ["rest","R1","R2"], ["rest","10kPa","1kPa"] };
    for g = 1:numel(row2Groups)
        nexttile; hold on; grid on;

        groupLevels = row2Groups{g};
        for lev = 1:numel(groupLevels)
            mask = string(G.diseaseCat)=="CT" & string(G.extrainfoCat)==groupLevels(lev);
            if any(mask)
                vals = G.Thickness{mask};
                if isempty(vals), continue; end
                prob = histcounts(double(vals), edges, 'Normalization','probability');
                prob(prob==0) = epsVal;
                plot(binCenters, prob, '-', 'Color', getColor(groupLevels(lev)), 'LineWidth', 1.8, 'DisplayName', groupLevels(lev));
            end
        end
        set(gca,'YScale','log','XLim',[0 xlimval]);
        xlabel('Thickness'); ylabel('Prevalence');
        title(sprintf('CT combined: %s', strjoin(groupLevels, '+')));
        legend('Location','best');
    end
    % Fill empty tiles in row2 (3,4,5) to prevent CCM falling in row2
    for k = 3:5
        nexttile; axis off;
    end

    % -------- Row3: CCM multi-extraLevel --------
    row3Groups = { ["rest","R1","R2"], ["rest","10kPa","1kPa"] };
    for g = 1:numel(row3Groups)
        nexttile; hold on; grid on;

        groupLevels = row3Groups{g};
        for lev = 1:numel(groupLevels)
            mask = string(G.diseaseCat)=="CCM" & string(G.extrainfoCat)==groupLevels(lev);
            if any(mask)
                vals = G.Thickness{mask};
                if isempty(vals), continue; end
                prob = histcounts(double(vals), edges, 'Normalization','probability');
                prob(prob==0) = epsVal;
                plot(binCenters, prob, '-', 'Color', getColor(groupLevels(lev)), 'LineWidth', 1.8, 'DisplayName', groupLevels(lev));
            end
        end
        set(gca,'YScale','log','XLim',[0 xlimval]);
        xlabel('Thickness'); ylabel('Prevalence');
        title(sprintf('CCM combined: %s', strjoin(groupLevels, '+')));
        legend('Location','best');
    end

    sgtitle('Thickness distributions (log Y-axis)');

end

% Color mapping as a function
function c = getColor(extraLevel)
    switch char(extraLevel)  % ensure char
        case 'rest',  c = [0 0 0];       % black
        case 'R1',    c = [0.6 0.7 1];   % light purple
        case 'R2',    c = [0.5 0.1 0.7]; % dark purple
        case '10kPa',  c = [0.6 1 0.6];  % light green
        case '1kPa',   c = [0 0.5 0];    % dark green
        otherwise,   c = [0.5 0.5 0.5]; % gray for safety
    end
end
