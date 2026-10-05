function T = Calc_Cell_AveragingMetrics(T, Input)
% Computes averaged CELL metrics (based on CellMask-derived measurements)
%
% REQUIRED INPUT COLUMNS:
%   AreaCell, MajAxCell, MinMajRatioCell, OrientationsCell
%
% INPUT:
%   Input.PixSize  - pixel size in µm
%
% OUTPUT columns added to T:
%   CellCount
%   CellAvPixSize
%   CellAvUmSize
%   CellAv_MajAx
%   CellAv_MajAx_um
%   CellAv_MinMajRatio
%   CellMeanAngle
%   CellCircVar
%   CellAlignIndex

PixSize_um = Input.PixSize;

% ── Validation ───────────────────────────────────────────────
requiredCols = ["AreaCell","MajAxCell","MinMajRatioCell","OrientationsCell"];

missing = requiredCols(~ismember(requiredCols, T.Properties.VariableNames));
if ~isempty(missing)
    error('Calc_Cell_AveragingMetrics: Missing columns: %s', strjoin(missing, ', '));
end

n = height(T);

% ── Preallocation ────────────────────────────────────────────
CellCount            = nan(n,1);
CellAvPixSize        = nan(n,1);
CellAvUmSize         = nan(n,1);
CellAv_MajAx         = nan(n,1);
CellAv_MajAx_um      = nan(n,1);
CellAv_MinMajRatio   = nan(n,1);
CellMeanAngle        = nan(n,1);
CellCircVar          = nan(n,1);
CellAlignIndex       = nan(n,1);

% ── Helper function ─────────────────────────────────────────
function v = getVec(tbl, col, row)
    raw = tbl.(col){row};
    if isnumeric(raw) && isvector(raw) && ~isempty(raw)
        v = raw(:);
    else
        v = [];
    end
end

% ── Main loop ───────────────────────────────────────────────
for i = 1:n

    vArea   = getVec(T,'AreaCell',i);
    if isempty(vArea), continue; end

    vMajAx  = getVec(T,'MajAxCell',i);
    vRatio  = getVec(T,'MinMajRatioCell',i);
    vAngle  = getVec(T,'OrientationsCell',i);

    % ── Count ───────────────────────────────────────────────
    CellCount(i) = sum(~isnan(vArea));

    if CellCount(i) >= 1

        % ── Size metrics ────────────────────────────────────
        CellAvPixSize(i) = mean(vArea,'omitnan');
        CellAvUmSize(i)  = CellAvPixSize(i) * PixSize_um^2;

        % ── Major axis ──────────────────────────────────────
        if ~isempty(vMajAx)
            CellAv_MajAx(i)    = mean(vMajAx,'omitnan');
            CellAv_MajAx_um(i) = CellAv_MajAx(i) * PixSize_um;
        end

        % ── Shape ratio ─────────────────────────────────────
        if ~isempty(vRatio)
            CellAv_MinMajRatio(i) = mean(vRatio,'omitnan');
        end

        % ── Orientation metrics (circular stats) ────────────
        valid = ~isnan(vAngle);

        if any(valid)
            theta = deg2rad(vAngle(valid)) * 2;

            C = mean(cos(theta));
            S = mean(sin(theta));

            R = sqrt(C^2 + S^2);   % alignment strength

            % mean angle (0–180°)
            meanAngle = rad2deg(atan2(S,C)) / 2;
            if meanAngle < 0
                meanAngle = meanAngle + 180;
            end

            CellMeanAngle(i)  = meanAngle;
            if CellCount(i) >= 3
                CellCircVar(i)    = 1 - R;
                CellAlignIndex(i) = R;
            else
                continue
            end
        end

    else
        continue
    end
end

% ── Attach to table ─────────────────────────────────────────
T.CellCount          = CellCount;
T.CellAvPixSize      = CellAvPixSize;
T.CellAvUmSize       = CellAvUmSize;
T.CellAv_MajAx       = CellAv_MajAx;
T.CellAv_MajAx_um    = CellAv_MajAx_um;
T.CellAv_MinMajRatio = CellAv_MinMajRatio;
T.CellMeanAngle      = CellMeanAngle;
T.CellCircVar        = CellCircVar;
T.CellAlignIndex     = CellAlignIndex;

end