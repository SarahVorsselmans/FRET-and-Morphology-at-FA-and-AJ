function T = Calc_AdhesionSkelCategoryStats(T, params)
% Calc_AdhesionSkelCategoryStats
% Per-category skeleton-pixel Int/FRET statistics.
%
% Uses the skeleton-pixel arrays already computed by Calc_Thickness
% (Skel_Int, Skel_FRET, skel_coords) and the category label image already
% computed by Calc_classifyAdhesions (MaskClass). For each skeleton pixel,
% looks up its category from MaskClass, then groups Skel_Int / Skel_FRET
% values by category — same logic as the skeleton-pixel fraction
% calculation in Calc_classifyAdhesions, but for intensity/FRET instead
% of area fraction.
%
% Requires T.MaskClass, T.skel_coords, T.Skel_Int (and T.Skel_FRET if FRET)
% to already exist.
%
% OUTPUT columns added to T:
%   SkelInt_Cat1, SkelInt_Cat2, ...    – cell, per-row vector of skeleton
%                                        Int values belonging to that category
%   SkelFRET_Cat1, SkelFRET_Cat2, ...  – same, for FRET (if params.FRET)
%   FA_summarySkelInt                  – row x nCats, mean Int per category
%   FA_summarySkelFRET                 – row x nCats, mean FRET per category (if FRET)

arguments
    T              table
    params.FRET    logical = false
    params.MIN_PX  (1,1) double = 30   % min skeleton pixels per category before summary computed
end

nRows = height(T);

% Derive nCats from MaskClass across all rows
nCats = 0;
for r = 1 : nRows
    mc = T.MaskClass{r};
    if ~isempty(mc)
        nCats = max(nCats, max(mc(:)));
    end
end

if nCats == 0
    warning('Calc_AdhesionSkelCategoryStats: no non-empty MaskClass found, nothing to do.');
    return
end

summarySkelInt  = nan(nRows, nCats);
summarySkelFRET = nan(nRows, nCats);
summarySkelThickness = nan(nRows, nCats);

for c = 1 : nCats
    T.("SkelInt_Cat" + c) = repmat({[]}, nRows, 1);
    T.("SkelThickness_Cat" + c) = repmat({[]}, nRows, 1);
    if params.FRET
        T.("SkelFRET_Cat" + c) = repmat({[]}, nRows, 1);
    end
end

allMaskClass  = T.MaskClass;
allSkelCoords = T.skel_coords;
allSkelInt    = T.Skel_Int;
allSkelThick   = T.Thickness_Skeleton; 
if params.FRET
    allSkelFRET = T.Skel_FRET;
end

fprintf('Calc_AdhesionSkelCategoryStats: processing %d rows...\n', nRows);

for row = 1 : nRows

    if mod(row, 50) == 0
        fprintf('  row %d / %d\n', row, nRows);
    end

    MaskClass = allMaskClass{row};
    coords    = allSkelCoords{row};
    skelInt   = allSkelInt{row};
    skelThick = allSkelThick{row};

    if isempty(MaskClass) || isempty(coords) || isempty(skelInt) || isempty(skelThick)
        continue
    end

    if params.FRET
        skelFRET = allSkelFRET{row};
        if isempty(skelFRET)
            continue
        end
    end

    nSkel = size(coords, 1);
    if numel(skelInt) ~= nSkel || numel(skelThick) ~= nSkel
        warning('Row %d: Skel_Int/Thickness_Skeleton length (%d) does not match skel_coords (%d), skipping.', ...
                 row, numel(skelInt), nSkel);
        continue
    end

    %% --- Look up category at each skeleton pixel coordinate ---
    sz       = size(MaskClass);
    linIdx   = sub2ind(sz, coords(:,1), coords(:,2));
    skelCat  = MaskClass(linIdx);   % nSkel x 1, category id per skeleton pixel (0 = background/unclassified)

    %% --- Group by category ---
    for c = 1 : nCats
        sel = (skelCat == c);
        if ~any(sel), continue; end

        T.("SkelInt_Cat" + c){row} = skelInt(sel);
        if numel(skelInt(sel)) >= params.MIN_PX
            summarySkelInt(row, c) = mean(skelInt(sel), 'omitnan');
        end

        T.("SkelThickness_Cat" + c){row} = skelThick(sel);          % <-- NEW
        if numel(skelThick(sel)) >= params.MIN_PX                   % <-- NEW
            summarySkelThickness(row, c) = mean(skelThick(sel), 'omitnan');  % <-- NEW
        end  

        if params.FRET
            T.("SkelFRET_Cat" + c){row} = skelFRET(sel);
            if numel(skelFRET(sel)) >= params.MIN_PX
                summarySkelFRET(row, c) = mean(skelFRET(sel), 'omitnan');
            end
        end
    end
end

T.FA_summarySkelInt = summarySkelInt;
T.FA_summarySkelThickness = summarySkelThickness;
if params.FRET
    T.FA_summarySkelFRET = summarySkelFRET;
end

fprintf('Done.\n');

end