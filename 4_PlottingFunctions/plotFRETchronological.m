function plotFRETchronological(T)
% PLOTFRETCHRONOLOGICAL Plot FRET% data chronologically for a selected DataFolder
%
% Input:
%   T - Table with columns: DataFolder, SampleFolder, EntryNumber, and FRET data column

%% --- Step 1: Let user select a DataFolder ---
dataFolders = unique(T.DataFolder);

fprintf('Available DataFolders:\n');
for i = 1:length(dataFolders)
    fprintf('  [%d] %s\n', i, dataFolders{i});
end

idx = input('Select DataFolder number: ');
selectedFolder = dataFolders{idx};

%% --- Step 2: Filter table for selected DataFolder ---
mask = strcmp(T.DataFolder, selectedFolder);
Tsub = T(mask, :);

%% --- Step 3: Get unique SampleFolders ---
sampleFolders = unique(Tsub.SampleFolder);
nSamples = length(sampleFolders);

%% --- Step 4: Define color palette (Prism-style, matches figure) ---
colors = [
    0.85, 0.33, 0.33;   % red
    0.93, 0.60, 0.55;   % light salmon
    0.95, 0.75, 0.65;   % peach
    0.20, 0.20, 0.20;   % black
    0.40, 0.40, 0.40;   % dark grey
    0.60, 0.60, 0.60;   % medium grey
    0.00, 0.55, 0.55;   % teal
    0.55, 0.75, 0.85;   % light blue
    0.20, 0.45, 0.70;   % blue
    0.85, 0.65, 0.20;   % gold
];
% Repeat colors if more samples than colors defined
colors = repmat(colors, ceil(nSamples / size(colors,1)), 1);

%% --- Step 5: Plot ---
figure('Color','white','Position',[100 100 560 380]);
hold on;

legendEntries = cell(nSamples, 1);

for s = 1:nSamples
    sf = sampleFolders{s};
    smask = strcmp(Tsub.SampleFolder, sf);
    Tsample = Tsub(smask, :);

    % Sort by EntryNumber for chronological x-axis
    Tsample = sortrows(Tsample, 'EntryNumber');

    x = Tsample.EntryNumber;
    y = Tsample.FRETav;

    plot(x, y, ...
        '-o', ...
        'Color',           colors(s,:), ...
        'MarkerFaceColor', colors(s,:), ...
        'MarkerEdgeColor', colors(s,:), ...
        'MarkerSize',       6, ...
        'LineWidth',        1.2);

    legendEntries{s} = sf;
end

%% --- Step 6: Formatting ---
xlabel('FRET Acquisition', 'FontSize', 12);
ylabel('FRET%',            'FontSize', 12);
title(selectedFolder,      'FontSize', 12, 'Interpreter', 'none');

xlim([min(Tsub.EntryNumber)-0.5, max(Tsub.EntryNumber)+0.5]);
ylim([15, 40]);                  % adjust to your data range

xticks(min(Tsub.EntryNumber) : max(Tsub.EntryNumber));

ax = gca;
ax.Box = 'off';
ax.TickDir = 'out';
ax.FontSize = 11;
ax.XColor = 'k';
ax.YColor = 'k';

legend(legendEntries, 'Location','eastoutside', ...
    'FontSize', 9, 'Box','off', 'Interpreter','none');

hold off;
end