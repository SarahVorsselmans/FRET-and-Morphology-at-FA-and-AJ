function T = LOOPcomputeCellSize(T)

hasMask = ismember("Mask", T.Properties.VariableNames);
if ~hasMask
    error('No mask source found');
end

nRows = height(T);

% Initialize CellMask as a cell column with NaNs
%T.CellSize         = nan(nRows,1);
T.CellMask = cell(height(T),1);
%T.CellMask(:) = {NaN};

% % Ensure the output vars exist with the correct types
% T.CellSize         = nan(height(T),1);     % mean per image, in µm²
% T.CellSizePerCell  = cell(height(T),1);    % 1xN vector (stored in a cell per row)
% T.NumCells         = zeros(height(T),1);   % diagnostic
% T.Status           = strings(height(T),1); % 'ok' or error message

for i = 1:nRows

    B = T.Mask{i};

    % Safety: ensure logical and non-empty
    if isa(B,'dip_image'), B = dip_array(B); end
    if isempty(B)
        %T.Status(i) = "empty_mask";
        continue
    end
    B = B ~= 0;

    % Call your existing function with figures OFF and I0=[]
    CellMask = computeCellSize(B);

    % % Write outputs
    % if isfield(out,'areas_um2') && ~isempty(out.areas_um2)
    %     T.CellSize(i)        = out.mean_um2;
    %     T.CellSizePerCell{i} = out.areas_um2(:)';  % 1xN row vector
    %     T.NumCells(i)        = numel(out.areas_um2);
    %     T.Status(i)          = "ok";
    % elseif isfield(out,'areas_px') && ~isempty(out.areas_px) && isempty(px)
    %     % If px was empty: store px²
    %     T.CellSize(i)        = mean(out.areas_px);
    %     T.CellSizePerCell{i} = out.areas_px(:)';
    %     T.NumCells(i)        = numel(out.areas_px);
    %     T.Status(i)          = "ok_no_px";
    % else
    %     T.Status(i) = "no_cells_detected";
    % end

    T.CellMask{i} = CellMask;

end

