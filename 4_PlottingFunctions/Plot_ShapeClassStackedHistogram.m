function Plot_ShapeClassStackedHistogram(T, options)
% -------------------------------------------------------------------------
% Generalized stacked histogram with 4 FA shape classes:
%   Cat1 | Cat2 | Cat3 | Cat4
%
% Adaptations:
%   • extrainfoCat levels resolved from table (stable order)
%   • optional user override via options.Levels
%   • flexible tiledlayout sizing
% -------------------------------------------------------------------------

%% ── ARGUMENTS ───────────────────────────────────────────────────────────
arguments
    T          table
    options.FRET       logical = false
    options.includePatterns cell = {}
    options.excludePatterns cell = {}
    options.Levels          cell = {}
    options.legendNames     cell = {}
end

FRET = options.FRET;
nCats = size(T.FA_summaryFrac, 2);  % calculates number of columns (= categories)

useIntCat = all(ismember("Int_Cat" + (1:nCats), T.Properties.VariableNames));

% rename
if isempty(options.legendNames)
    options.legendNames = {'ThinLine','EllipsoidLine','SmallRound','BigRound'};
end

%% ── FILTER ROWS ─────────────────────────────────────────────────────────
T = Filter_IncludeExclude(T, options.includePatterns, options.excludePatterns);

%% ── RESOLVE DISEASE + EXTRAINFO LEVELS ──────────────────────────────────
diseaseList = ["CT","CCM"];
nDisease    = numel(diseaseList);

% Resolve extrainfoCat levels from table (stable order)
if iscategorical(T.extrainfoCat)
    allLevels = cellstr(unique(T.extrainfoCat, 'stable'));
else
    allLevels = unique(T.extrainfoCat, 'stable');
    if isstring(allLevels)
        allLevels = cellstr(allLevels);
    end
end

% Apply optional user override
if ~isempty(options.Levels)
    bad = options.Levels(~ismember(options.Levels, allLevels));
    if ~isempty(bad)
        error('Requested extrainfoCat level(s) not found: %s', ...
              strjoin(bad, ', '));
    end
    infoList = options.Levels;
else
    infoList = allLevels;
end

nInfo = numel(infoList);

%% ── FIGURE & LAYOUT ─────────────────────────────────────────────────────
figure;
tiledlayout(nDisease, nInfo, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% ── MAIN PANEL LOOP ─────────────────────────────────────────────────────
for d = 1:nDisease
    for e = 1:nInfo

        nexttile;
        hold on;

        mask = T.diseaseCat == diseaseList(d) & ...
               T.extrainfoCat == infoList{e};
        subT = T(mask, :);

        % Counters
        % nCat1 = 0; nCat2 = 0; nCat3 = 0; nCat4 = 0;

        % FRET accumulation
        % FRET_Cat1 = []; FRET_Cat2 = []; FRET_Cat3 = []; FRET_Cat4 = [];
        nPerCat = zeros(1, nCats);
        FRET_perCat = cell(1, nCats);
        for c = 1:nCats
            FRET_perCat{c} = [];
        end
        
        if useIntCat
            % Count localizations from Int_Cat
            for i = 1:height(subT)
                for c = 1:nCats
                    colInt = "Int_Cat" + c;
                    vals = subT.(colInt){i};
                    nPerCat(c) = nPerCat(c) + numel(vals);
                    if FRET
                        colFRET = "FRET_Cat" + c;
                        FRET_perCat{c} = [FRET_perCat{c}; subT.(colFRET){i}];
                    end
                end
            end
        else
            % Sum fractions from FA_summaryFrac
            if height(subT) > 0
                nPerCat = sum(subT.FA_summaryFrac, 1, 'omitnan');  % 1 x nCats
            end
            if FRET
                for i = 1:height(subT)
                    for c = 1:nCats
                        colFRET = "FRET_Cat" + c;
                        FRET_perCat{c} = [FRET_perCat{c}; subT.(colFRET){i}];
                    end
                end
            end
        end

        %% ── STATS ────────────────────────────────────────────────────────
        total = sum(nPerCat);
        if total > 0
            perc = nPerCat / total * 100;
        else
            perc = zeros(1, nCats);
        end
        
        if FRET
            FRETmean = cellfun(@(v) mean(v, 'omitnan'), FRET_perCat);
        end

        %% ── PLOT ─────────────────────────────────────────────────────────
        h = bar(1, perc, 'stacked');
        for c = 1:nCats
            h(c).FaceColor = getColorClass(c);
        end

        ylim([0 100]);
        xlim([0.5 1.5]);
        set(gca,'XTick',[]);
        % grid on;

        if e == 1
            ylabel('Fraction (%)');
        end

        % Annotations
        yBottom = 0;
        for c = 1:nCats
            if perc(c) > 0
                yCenter = yBottom + perc(c)/2;
                if FRET
                    txt = sprintf('%.1f%%  FRET=%.2f', perc(c), FRETmean(c));
                else
                    txt = sprintf('%.1f%%', perc(c));
                end
                text(1, yCenter, txt, ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', ...
                    'Color','w', 'FontWeight','bold');
            end
            yBottom = yBottom + perc(c);
        end

        nTotal = total;

        if useIntCat
            title(sprintf('%s – %s\n(FAs = %d,  N = %d)', ...
                diseaseList(d), infoList{e}, nTotal, height(subT)), ...
                'FontSize',10,'FontWeight','bold');
        else
            title(sprintf('%s – %s\n(N = %d)', ...
                diseaseList(d), infoList{e}, height(subT)), ...
                'FontSize',10,'FontWeight','bold');
        end
    end
end

%% ── LEGEND ──────────────────────────────────────────────────────────────
if isempty(options.legendNames)
    options.legendNames = arrayfun(@(c) "Cat"+c, 1:nCats, 'UniformOutput', false);
end
legend(options.legendNames, 'Location', 'bestoutside');
end