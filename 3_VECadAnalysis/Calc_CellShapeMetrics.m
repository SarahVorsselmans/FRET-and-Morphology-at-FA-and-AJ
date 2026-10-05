function T = Calc_CellShapeMetrics(T)
% Calculate shape metrics based on T.CellMask instead of T.Mask
%
% REQUIRED INPUT:
%   T.CellMask{i}  - binary mask of cell(s)
%
% OUTPUT columns added to T:
%   .AreaCell
%   .MajAxCell
%   .MinAxCell
%   .MinMajRatioCell
%   .OrientationsCell

% ── Validation ────────────────────────────────────────────────────────────
if ~ismember("CellMask", T.Properties.VariableNames)
    error('Calc_CellShapeMetrics: T.CellMask not found.');
end

nRows = height(T);

% ── Preallocate output columns ────────────────────────────────────────────
T.AreaCell            = cell(nRows,1);
T.MajAxCell           = cell(nRows,1);
T.MinAxCell           = cell(nRows,1);
T.MinMajRatioCell     = cell(nRows,1);
T.OrientationsCell    = cell(nRows,1);

% ── Main loop ─────────────────────────────────────────────────────────────
for i = 1:nRows

    BW = T.CellMask{i};

    % Ensure logical
    BW = logical(BW);

    % Connected components (cells)
    CC   = bwconncomp(BW, 4);
    nObj = CC.NumObjects;

    if nObj == 0
        T.AreaCell{i}         = [];
        T.MajAxCell{i}        = [];
        T.MinAxCell{i}        = [];
        T.MinMajRatioCell{i}  = [];
        T.OrientationsCell{i} = [];
        continue
    end

    % ── regionprops ──────────────────────────────────────────────
    stats = regionprops(CC, ...
        'Area', ...
        'MajorAxisLength', ...
        'MinorAxisLength', ...
        'Orientation');

    % ── Extract metrics ──────────────────────────────────────────
    Area  = vertcat(stats.Area);
    MajAx = vertcat(stats.MajorAxisLength);
    MinAx = vertcat(stats.MinorAxisLength);

    MinMajRatio = MinAx ./ MajAx;

    % Orientation filtering (same logic as original)
    aspectRatio = MajAx ./ max(MinAx, 1);
    keep        = aspectRatio > 1.5;

    angles_deg  = nan(nObj,1);
    rawAngles   = vertcat(stats.Orientation);
    angles_deg(keep) = rawAngles(keep);

    % ── Store results ────────────────────────────────────────────
    T.AreaCell{i}         = Area;
    T.MajAxCell{i}        = MajAx;
    T.MinAxCell{i}        = MinAx;
    T.MinMajRatioCell{i}  = MinMajRatio;
    T.OrientationsCell{i} = angles_deg;

end
end
