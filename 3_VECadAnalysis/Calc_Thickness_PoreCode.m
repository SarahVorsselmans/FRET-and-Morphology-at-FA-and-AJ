function T = computeThickness_PoreCode(T)
% computeThickness_PoreCode: Calculates junction thickness and plots distribution
% Inputs:
%   T.FRETimageFA/FRETindex - (converted to) logical 2D array (segmentation mask)
% Outputs:
%   T.ThicknessBub   - per thickness value
%   T.ThicknessBubAv - mean thickness (µm/pixels??)

    useMask = ismember("FRETimageFA", T.Properties.VariableNames);
    useFRET = ismember("FRETindex", T.Properties.VariableNames);
    
    if ~useMask && ~useFRET
        warning('Neither Mask nor FRETindex exist. Nothing to compute.');
        return
    end
    
    % Parameters
    % pxSize      = 0.24; % µm per pixel (adjust if needed)
    % binWidth_um = 0.5; % bin width for histogram-line plot (µm)
    
    % Preallocate: one cell per row for the diameter vectors
    T.ThicknessBub   = cell(height(T), 1);
    % Preallocate: one scalar per row for the averages
    T.ThicknessBubAv = nan(height(T), 1);

    for i = 1:height(T) % = nRows
        % Choose source of mask
        if useMask && ~isempty(T.FRETimageFA{i})
            % User-supplied binary mask (0/1)
            BW = ~isnan(T.FRETimageFA{i});
            % BW = logical(BW);
        elseif useFRET && ~isempty(T.FRETindex{i})
            % Guard against missing/empty masks
            BW = T.FRETindex{i};
            % Replace invalid values with background (100)
            BW(~isfinite(BW)) = 100;
            % Build binary mask: 100 -> 0 (background), anything else -> 1 (foreground)
            BW = BW ~= 100;
        end
        if isempty(BW)
            T.ThicknessBub{i}   = [];      % keep empty
            T.ThicknessBubAv(i) = NaN;     % average is NaN
            continue
        end

        % Ensure logical for bwmorph (it expects a binary/logical image)
        BW = logical(BW);

        %%% NEEED TO INSERT: BETTER FILLING OF HOLES - SEE AREA AND SKELETON CODE

        % Run your function on the closed mask
        pores2D = getPoreProps2D(bwmorph(BW, 'close'));

        % Get the diameter vector (Nx1 double). Make sure it's a column.
        diam = pores2D.diameter(:);

        % Store the vector into a cell and the average into a numeric column
        T.ThicknessBub{i}   = diam;
        T.ThicknessBubAv(i) = mean(diam, 'omitnan');
    end
end