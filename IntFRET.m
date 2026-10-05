function T = IntFRET(data)
%COMPUTEFRET Compute FRET ratio image and summary statistics.
%
% INPUTS
%   mask : 2D binary/uint/logical array (background = 0). Objects are the non-zero regions.
%   DD  : 2D array (same size as mask). Pixel values typically 0–255.
%   DA  : 2D array (same size as mask). Pixel values typically 0–255.
%
% CALCULATIONS
%   1) Segment DD and DA by the mask.
%   2) Compute FRET image as DA ./ DD (element-wise) only inside mask.
%      - Outside mask -> NaN.
%      - Where DD == 0 -> NaN to avoid Inf/Divide-by-zero.
%   3) Find connected components (separate objects) in mask and compute
%      the average FRET per object (omit NaNs).
%   4) Compute the average FRET for the whole image within the mask.
%
% OUTPUT
%   T : 1x1 table with columns:
%       - 'mask'        : 1x1 cell containing the mask array
%       - 'DD'         : 1x1 cell containing the DD array
%       - 'DA'         : 1x1 cell containing the DA array
%       - 'FRETimage'   : 1x1 cell containing the FRET ratio image (double, NaN outside mask)
%       - 'FRETperObject' : 1x1 cell containing an Nx1 double vector (mean FRET per object)
%       - 'FRETav'      : scalar double (mean FRET across all masked pixels)
%
% NOTES
%   - This function does not assume a fixed size (e.g., 512x512); it only
%     requires that all inputs have the same size.
%   - Requires Image Processing Toolbox for bwconncomp (preferred). If not
%     available, you can replace with bwlabel.
%
% EXAMPLE
%   T = computeFRET(mask, DD, DA);
%   imagesc(T.FRETimage{1}); axis image; colorbar; title('FRET ratio');
%
% Author: <your name>, <date>

    % -------------------------
    % Inputs
    % -------------------------
    % Coerce types
    mask = data.mask ~= 0;              % logical mask (non-zero treated as true)
    DD  = data.DD;
    DA  = data.DA;

    % --- Median filter channels with a 9x9 window ---
    win = [9 9];
    DD_filt = medfilt2(DD, win, 'symmetric');
    DA_filt = medfilt2(DA, win, 'symmetric');

    % -------------------------
    % Compute FRET image (DA ./ DD) inside mask
    % -------------------------
    % FRETimage = nan(sz);           % NaN outside mask by default
    % denom = DD;
    % denom(denom == 0) = NaN;       % avoid division by zero (-> NaN)
    % % Compute ratio only where mask is true
    % idx = mask & ~isnan(denom) & ~isnan(DA);
    % FRETimage(idx) = DA(idx) ./ denom(idx);

    FRETimage = nan(size(mask));
    denom = DD_filt;
    denom(denom == 0) = NaN;
    
    idx = mask & ~isnan(denom) & ~isnan(DA_filt);
    FRETimage(idx) = DA_filt(idx) ./ denom(idx);

    % -------------------------
    % Per-object averages (connected components)
    % -------------------------
    % 8-connectivity is common for cell-like objects; change to 4 if preferred
    try
        CC = bwconncomp(mask, 8);
        pixelLists = CC.PixelIdxList;
    catch
        % Fallback if bwconncomp not available (older MATLAB)
        [L, num] = bwlabel(mask, 8);
        pixelLists = cell(num,1);
        for k = 1:num
            pixelLists{k} = find(L == k);
        end
        CC.NumObjects = num;
    end

    nObj = numel(pixelLists);
    FRETperObject = nan(nObj, 1);
    for k = 1:nObj
        vals = FRETimage(pixelLists{k});
        FRETperObject(k) = mean(vals, 'omitnan');
    end

    % -------------------------
    % Whole-image (masked) average
    % -------------------------
    FRETav = mean(FRETimage(mask), 'omitnan');

    % -------------------------
    % Package into a 1-row table
    % -------------------------
    % Use cells to store matrices/vectors so this table can be vertically
    % concatenated with other images later on.
    T = table( ...
        {mask}, {DD}, {DA}, {FRETimage}, {FRETperObject}, FRETav, ...
        'VariableNames', {'mask','DD','DA','FRETimage','FRETperObject','FRETav'} ...
    );
end