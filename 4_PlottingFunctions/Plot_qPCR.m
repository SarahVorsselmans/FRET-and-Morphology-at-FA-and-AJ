%% qPCR Analysis - CCM2, ROCK1, ROCK2 (Vin & VE combined per condition)
% Converted from Python. Changes vs. the original:
%   - Only CCM2 / ROCK1 / ROCK2 are plotted (all 6 targets are still used
%     for the DeltaCt/DeltaDeltaCt calculation, since VE CT normalizes
%     every target).
%   - Vin and VE data for the same condition (e.g. "Vin CT" + "VE CT")
%     are pooled into a single bar (mean +/- SD across both cell types).
%   - Individual Vin/VE data points are overlaid as a scatter on top of
%     the bars: circles = Vin, triangles = VE.
%   - Bars are wider (BarWidth below).
%   - Coloring (blue = CT, red = CCM, darker = fewer replicates) is kept.

close all; clear; clc;

%% ---------------------- 1. LOAD DATA ----------------------
dataPathRep = 'C:\Users\sarah\OneDrive - KU Leuven\Results_Grenoble\qPCR\INVADE _20251202_185547_796BR03005_sarahv_values_resaved.xlsx';
dfRep = readtable(dataPathRep);

% Keep only the relevant columns
dfRep = dfRep(:, {'Target', 'Sample', 'Cq'});

% Standardize text columns to string (readtable may import as cellstr)
dfRep.Target = string(dfRep.Target);
dfRep.Sample = string(dfRep.Sample);

%% ---------------------- 2. REFERENCE GENES ----------------------
ref_targets = ["GAPDH", "RANBP1"];
dfRef = dfRep(ismember(dfRep.Target, ref_targets), :);

% Average Cq per Sample + Target
dfRef_avg = groupsummary(dfRef, {'Sample', 'Target'}, 'mean', 'Cq');
dfRef_avg.Properties.VariableNames{'mean_Cq'} = 'Cq_ref_avg';

% Average across reference targets -> one value per sample
dfRef_final = groupsummary(dfRef_avg, 'Sample', 'mean', 'Cq_ref_avg');
dfRef_final.Properties.VariableNames{'mean_Cq_ref_avg'} = 'Cq_REF';
dfRef_final = dfRef_final(:, {'Sample', 'Cq_REF'});

disp(dfRef_final)

%% ---------------------- 3. DELTA CT ----------------------
target_genes = ["CCM2", "ROCK1", "ROCK2", "Vinculin", "VeCad", "HPRT1"];
dfTargets = dfRep(ismember(dfRep.Target, target_genes), :);

dfTargets = outerjoin(dfTargets, dfRef_final, ...
    'Keys', 'Sample', 'MergeKeys', true, 'Type', 'left');

dfTargets.DeltaCt = dfTargets.Cq - dfTargets.Cq_REF;

%% ---------------------- 4. DELTA DELTA CT & EXPRESSION ----------------------
control_sample = "VE CT";

ctrlRows = dfTargets(dfTargets.Sample == control_sample, :);
control_means = groupsummary(ctrlRows, 'Target', 'mean', 'DeltaCt');
control_means.Properties.VariableNames{'mean_DeltaCt'} = 'DeltaCt_control';
control_means = control_means(:, {'Target', 'DeltaCt_control'});

dfTargets = outerjoin(dfTargets, control_means, ...
    'Keys', 'Target', 'MergeKeys', true, 'Type', 'left');

dfTargets.DeltaDeltaCt = dfTargets.DeltaCt - dfTargets.DeltaCt_control;
dfTargets.Expression   = 2 .^ (-dfTargets.DeltaDeltaCt);

%% ---------------------- 5. SPLIT SAMPLE INTO CellType + Condition ----------------------
% "Vin CT R1" -> CellType = "Vin", Condition = "CT R1"
% "VE CCM"    -> CellType = "VE",  Condition = "CCM"
CellType  = strings(height(dfTargets), 1);
Condition = strings(height(dfTargets), 1);

for i = 1:height(dfTargets)
    s = dfTargets.Sample(i);
    if startsWith(s, "Vin ")
        CellType(i)  = "Vin";
        Condition(i) = extractAfter(s, "Vin ");
    elseif startsWith(s, "VE ")
        CellType(i)  = "VE";
        Condition(i) = extractAfter(s, "VE ");
    else
        CellType(i)  = "";
        Condition(i) = "";
    end
end

dfTargets.CellType  = CellType;
dfTargets.Condition = Condition;

%% ---------------------- 6. FILTER TO GENES / CONDITIONS OF INTEREST ----------------------
OrderTargetPlot = ["CCM2", "ROCK1", "ROCK2"];
OrderCondition  = ["CT", "CT R1", "CT R2", "CCM", "CCM R1", "CCM R2"];

dfPlot = dfTargets(ismember(dfTargets.Target, OrderTargetPlot) & ...
                    ismember(dfTargets.Condition, OrderCondition), :);

%% ---------------------- 7. BUILD BAR + SCATTER DATA WITH SPACING ----------------------
barX      = [];
barHeight = [];
barSD     = [];
barColor  = [];
barLabel  = strings(0, 1);

targetLabelPos  = [];
targetLabelName = strings(0, 1);
targetGroupStart = [];
targetGroupEnd   = [];

scatterX     = [];
scatterY     = [];
scatterShape = strings(0, 1);

intraGap = 1;   % gap between the CT sub-block and CCM sub-block, within a target
interGap = 2;   % gap between targets
jitterW  = 0.12;

pos = 0;
for t = 1:numel(OrderTargetPlot)
    tgt = OrderTargetPlot(t);
    groupStart = pos;

    for c = 1:numel(OrderCondition)
        cond = OrderCondition(c);

        rows = dfPlot.Target == tgt & dfPlot.Condition == cond;
        vals = dfPlot.Expression(rows);
        cellTypes = dfPlot.CellType(rows);

        if isempty(vals)
            continue
        end

        barX(end+1)       = pos;                       %#ok<*SAGROW>
        barHeight(end+1)  = mean(vals, 'omitnan');
        barSD(end+1)      = std(vals, 'omitnan');
        barColor(end+1,:) = getColor(cond);
        barLabel(end+1)   = cond;

        for k = 1:numel(vals)
            if cellTypes(k) == "Vin"
                baseOffset = -0.18;
                shape = "o";
            else
                baseOffset = 0.18;
                shape = "^";
            end
            scatterX(end+1)     = pos + baseOffset + (rand - 0.5) * jitterW;
            scatterY(end+1)     = vals(k);
            scatterShape(end+1) = shape;
        end

        pos = pos + 1;

        % gap between the CT and CCM sub-blocks within this target
        if cond == "CT R2"
            pos = pos + intraGap;
        end
    end

    groupEnd = pos - 1;
    targetGroupStart(end+1) = groupStart;
    targetGroupEnd(end+1)   = groupEnd;
    targetLabelPos(end+1)   = (groupStart + groupEnd) / 2;
    targetLabelName(end+1)  = tgt;

    pos = pos + interGap;
end

%% ---------------------- 8. PLOT ----------------------
barWidth     = 0.85;                          % wider bars
yMaxEstimate = max(barHeight + barSD) * 1.3;

figure('Position', [100 100 1300 600]);
hold on

% --- Shaded background per target group ---
for t = 1:numel(OrderTargetPlot)
    xStart = targetGroupStart(t) - 0.5;
    xEnd   = targetGroupEnd(t)   + 0.5;
    patch([xStart xEnd xEnd xStart], [0 0 yMaxEstimate yMaxEstimate], ...
          [0.85 0.85 0.85], 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
          'HandleVisibility', 'off');
end

% --- Bars (plotted one at a time so each can have its own alpha) ---
for i = 1:numel(barX)
    bar(barX(i), barHeight(i), barWidth, ...
        'FaceColor', barColor(i, 1:3), 'FaceAlpha', barColor(i, 4), ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% --- Error bars ---
errorbar(barX, barHeight, barSD, 'k', 'LineStyle', 'none', ...
         'CapSize', 5, 'HandleVisibility', 'off');

% --- Scatter overlay: Vin = circle, VE = triangle ---
% isVin = scatterShape == "o";
% isVE  = scatterShape == "^";
% 
% hVin = scatter(scatterX(isVin), scatterY(isVin), 45, 'k', 'o', ...
%                'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'w');
% hVE  = scatter(scatterX(isVE), scatterY(isVE), 45, 'k', '^', ...
%                'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'w');
% 
% legend([hVin hVE], {'Vinculin', 'VE-Cadherin'}, ...
%        'Location', 'northoutside', 'Orientation', 'horizontal');

% --- X ticks with condition names ---
xticks(barX);
xticklabels(barLabel);
xtickangle(45);

% --- Target labels above each group ---
for t = 1:numel(OrderTargetPlot)
    text(targetLabelPos(t), yMaxEstimate * 0.97, targetLabelName(t), ...
         'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
end

ylabel('Expression (mean \pm SD)');
title('qPCR Expression - CCM2 / ROCK1 / ROCK2 (Vin + VE combined)');
ylim([0 yMaxEstimate]);
xlim([-1, max(barX) + 1]);
box off
hold off

%% ---------------------- LOCAL FUNCTIONS ----------------------
function c = getColor(cond)
% Same coloring scheme as the original: blue = CT, red = CCM,
% darker/more opaque = the "main" replicate, lighter = R1/R2.
    if contains(cond, "CT R2")
        c = [0 0 1 0.3];
    elseif contains(cond, "CT R1")
        c = [0 0 1 0.6];
    elseif contains(cond, "CT")
        c = [0 0 1 1.0];
    elseif contains(cond, "CCM R2")
        c = [1 0 0 0.3];
    elseif contains(cond, "CCM R1")
        c = [1 0 0 0.6];
    elseif contains(cond, "CCM")
        c = [1 0 0 1.0];
    else
        c = [0.5 0.5 0.5 1.0];
    end
end