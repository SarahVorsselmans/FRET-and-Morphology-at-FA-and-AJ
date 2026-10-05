function T = Calc_classifyAdhesions(T, params)
% Calc_classifyAdhesions
% Classify adhesion segments into up to N categories defined externally.
%
% Categories are evaluated IN ORDER – once a fragment is assigned to a
% category it is excluded from all subsequent category tests.
%
% Built-in category names: Cat1, Cat2, Cat3 (extensible via AJCatCategories)
%   Cat1 (default green)  – smooth thin adhesions
%   Cat2 (default red)    – thick cheese-like adhesions
%   Cat3 (default blue)   – all remaining fragments
%
% ── External category definition ────────────────────────────────────────
% Define Input.AJCatCategories as a cell array OUTSIDE this function and
% pass it in via params.AJCatCategories.  Each element is a cell row that
% lists conditions as triplets {feature, operator, threshold} joined by
% implicit AND.  Example:
%
%   Input.AJCatCategories = {
%       { 'Thickness','<',8,  'StdRatio','<',0.35 },  ... % Cat1
%       { 'Thickness','>',16, 'Solidity','<',0.8  },  ... % Cat2
%   };                                                     % Cat3 = remainder
%
% Supported feature names (must match statsTable column names exactly):
%   SkelLength  Thickness  ThicknessStd  StdRatio
%   Tortuosity  BranchCount  EulerNum  Solidity  AspectRatio
%
% Supported operators:  '<'  '<='  '>'  '>='  '=='  '~='
%
% The LAST category in AJCatCategories is treated as a catch-all (like an
% else-branch) – it receives every fragment not yet assigned.  If only one
% entry is provided it acts as Cat1 and everything else becomes Cat2.
%
% ── Usage ────────────────────────────────────────────────────────────────
%   % Minimal – use built-in defaults
%   [lm, tbl, rgb] = Calc_classifyAdhesions(BW);
%
%   % Custom categories loaded from external cell
%   Input.AJCatCategories = {
%       { 'Thickness','<',8, 'StdRatio','<',0.35 },  ... % Cat1
%       { 'Thickness','>',16,'Solidity','<',0.8  },  ... % Cat2
%   };
%   [lm, tbl, rgb] = Calc_classifyAdhesions(BW, Input);
%
% INPUT
%   BW      : binary image (logical / double, non-zero = foreground)
%   params  : name-value or struct – all fields optional (see arguments block)
%
% OUTPUT
%   labelMap   : uint8 image  (0 = background, 1 = Cat1, 2 = Cat2, …)
%   statsTable : table with one row per fragment + computed features
%   RGB        : colour visualisation (Cat1=green, Cat2=red, Cat3=blue, …)

%% ── REMARKS! ────────────────────────────────────────────────────────────
% I deleted a lot of parameters to be calculated (e.g. SkelLength, Solidity
% etc.[see section 2]) to improve speed. If you want to use those uncomment
% them.

%% ── Default parameters ──────────────────────────────────────────────────
arguments
    %BW                              logical
    T                               table

    % ── Preprocessing ────────────────────────────────────────────────────
    params.min_area                 (1,1) double  = 50    % minimum CC area (px) before splitting – removes noise
    params.ws_h_minima              (1,1) double  = 1.0   % 1.5 % imhmin depth for watershed; larger → fewer splits
    params.min_frag_area            (1,1) double  = 30    % discard fragments smaller than this after splitting

    % ── Unused splitting params kept for future use ───────────────────────
    params.split_at_branches        logical       = true  % (reserved) cut skeleton at branch-points
    params.split_at_widthchange     logical       = true  % (reserved) cut at width-change transitions
    params.width_jump_ratio         (1,1) double  = 2.0   % (reserved) ratio threshold for width-change cuts
    params.width_smooth_sigma       (1,1) double  = 3     % (reserved) Gaussian σ for smoothing width profile
    params.min_split_area           (1,1) double  = 30    % (reserved) min area after skeleton-based split

    % ── Classification thresholds (used in built-in fallback categories) ──
    params.L_thr                    (1,1) double  = 40    % min skeleton length for "long" fragments
    params.T_low                    (1,1) double  = 5.0   % max thickness (px diameter) for Cat1 thin
    params.T_high                   (1,1) double  = 16    % min thickness (px diameter) for Cat2 thick
    params.tort_thr                 (1,1) double  = 1.25  % max tortuosity  (1.0 = perfectly straight)
    params.branch_thr               (1,1) double  = 2     % max branch-point count for Cat1
    params.solidity_thr             (1,1) double  = 0.75  % solidity below this → Cat2 (holey / lumpy)
    params.std_ratio_thr            (1,1) double  = 0.3   % max thickness StdDev/Mean for Cat1 (uniform width)

    % ── External category definitions ─────────────────────────────────────
    % Cell array – each element defines one category (see header above).
    % Leave empty {} to use the built-in Cat1/Cat2/Cat3 logic.
    params.AJCatCategories          cell          = {}
end

% % Remove DIPImage interference
% dipPath = fileparts(which('regionprops'));
% if contains(dipPath, 'DIPimage', 'IgnoreCase', true)
%     rmpath(dipPath);
% end
% % restore DIPImage location at end of function
% if ~isempty(dipPath)
%     addpath(dipPath);
% end

nRows = height(T);
nCats = numel(params.AJCatCategories); % Derive number of categories
summaryFrac = nan(nRows, nCats); % pre-allocate column
maskClasses = cell(nRows, 1);    % pre-allocate column

fprintf('Calc_classifyAdhesions: processing %d rows...\n', nRows);

allMasks = T.Mask; % Pull masks out of table once – avoids repeated table indexing in parfor

for row = 1 : nRows

    if mod(row, 50) == 0
        fprintf('  row %d / %d\n', row, nRows);
    end

    %% ── Pre-processing ───────────────────────────────────────────────────────
    BW = logical(allMasks{row});
    if ~any(BW(:)) % Guard: skip empty masks
        maskClasses{row} = zeros(size(BW));
        continue
    end

    BW = bwareaopen(BW, params.min_area);
    if ~any(BW(:)) % Guard: skip empty masks
        maskClasses{row} = zeros(size(BW));
        continue
    end

    BW = imclose(BW, strel('disk', 1)); % 2
    
    %% ── Distance transform ───────────────────────────────────────────────────
    % D = bwdist(~BW);
    
    %% ════════════════════════════════════════════════════════════════════════
    %  Watershed on the negated distance transform places cut-lines at saddle
    %  points (thin necks between thick blobs), which run PERPENDICULAR to
    %  the adhesion axis – exactly where morphological category changes occur.
    %
    %  ws_h_minima suppresses shallow minima in −D before the watershed so
    %  that near-flat regions are not over-segmented.
    %
    %  New: Watershed sectioning with hole preservation
    %  Holes (enclosed background regions) are identified and protected
    %  before watershed, then restored into each section afterward.
    % ════════════════════════════════════════════════════════════════════════
    % ── Step 1: Identify holes ───────────────────────────────────────────────
    % Holes = background pixels NOT connected to the image border.
    % imfill finds them; subtracting BW isolates only the enclosed background.
    BW_filled_all = imfill(BW, 'holes');
    holes_all     = BW_filled_all & ~BW;          % every enclosed background region

    % ── Step 1b: Keep only SMALL holes (likely gaps/artifacts, not real cells) ─
    max_hole_area = 500;   % max size of holes
    holes         = bwareaopen(holes_all, 1) & ~bwareaopen(holes_all, max_hole_area);
    % (bwareaopen(holes_all,1) just cleans noise; the second term removes
    %  components that survive an area-opening at max_hole_area, i.e. big ones)
    BW_filled = BW | holes;   % only re-fill the small holes, leave big ones as background
    
    % ── Step 2: Run watershed on a HOLE-FILLED foreground ───────────────────
    % Filling holes first means the distance transform sees solid blobs,
    % so saddle-point cuts land on thin necks – never inside a hole.
    D_filled   = bwdist(~BW_filled);      % distance transform of filled shape
    D_neg      = -D_filled;
    D_neg(~BW_filled) = 0;
    
    D_sup = imhmin(D_neg, params.ws_h_minima);
    L_ws  = watershed(D_sup);
    
    % ── Step 3: Map watershed labels back to original (unfilled) foreground ──
    % Watershed ran on BW_filled; restrict labels to the real foreground BW.
    BW_work = BW & (L_ws > 0);
    
    % ── Step 4: Re-punch holes into every section ────────────────────────────
    % Any foreground pixel that belongs to a hole is cleared so that
    % each section retains the hole pixels that fall within its extent.
    BW_work = BW_work & ~holes;
    
    % ── Step 5: Small-fragment removal (unchanged) ───────────────────────────
    BW_work = bwareaopen(BW_work, params.min_frag_area);
    if ~any(BW_work(:)) % Guard: skip if no fragments survive
        maskClasses{row} = zeros(size(BW));
        continue
    end
    
    %% ════════════════════════════════════════════════════════════════════════
    %  STEP 2  –  Feature extraction on split fragments
    % ════════════════════════════════════════════════════════════════════════
    % Fill holes in each connected component for skeleton-based metrics.
    % Solidity and EulerNumber are computed on the original (holey) shapes.
    BW_work_filled = imfill(BW_work, 'holes');
    
    Skel2  = bwskel(BW_work_filled);
    %BP2    = bwmorph(Skel2, 'branchpoints');
    D2     = bwdist(~BW_work_filled);
    
    CC     = bwconncomp(BW_work_filled, 8);
    % props  = regionprops(CC, 'Area','Perimeter','Eccentricity', ...
    %                          'MajorAxisLength','MinorAxisLength');
    
    % Solidity and EulerNumber are computed on the original shapes with holes.
    CC_holey    = bwconncomp(BW_work, 8);
    props_holey = regionprops(CC_holey, 'Solidity', 'EulerNumber');
    
    numObj = CC.NumObjects;
    
    % Pre-allocate feature vectors
    % length_skel   = zeros(numObj,1);
    thickness     = zeros(numObj,1);
    % thickness_std = zeros(numObj,1);
    % tortuosity    = nan(numObj,1);
    % branch_count  = zeros(numObj,1);
    % solidity_vec  = zeros(numObj,1);
    euler_vec     = zeros(numObj,1);
    % aspect_ratio  = zeros(numObj,1);
    % std_ratio     = zeros(numObj,1);
    
    for i = 1 : numObj
        mask = false(size(BW_work_filled));
        mask(CC.PixelIdxList{i}) = true;
    
        % sk = Skel2 & mask;
    
        % % Skeleton length (px)
        % length_skel(i) = sum(sk(:));
    
        % Thickness: diameter = 2 × mean of distance-transform radii
        radii = D2(mask);
        thickness(i)     = 2 * mean(radii);
        % thickness_std(i) = std(radii);
        % std_ratio(i)     = thickness_std(i) / (thickness(i) + eps);
    
        % % Tortuosity: skeleton arc-length / straight-line endpoint distance
        % ep = bwmorph(sk, 'endpoints');
        % [ey, ex] = find(ep);
        % if numel(ex) >= 2
        %     euclid        = hypot(ex(1)-ex(end), ey(1)-ey(end));
        %     tortuosity(i) = length_skel(i) / (euclid + eps);
        % end
    
        % % Number of skeleton branch-points within this fragment
        % branch_count(i) = sum(BP2(mask));
    
        % Solidity and EulerNumber: from original shapes that retain holes
        euler_vec(i)    = props_holey(i).EulerNumber;
        % solidity_vec(i) = props_holey(i).Solidity;
    
        % aspect_ratio(i) = props(i).MajorAxisLength / (props(i).MinorAxisLength + eps);
    end
    
    %% ── Assemble feature table (used both for classification and output) ─────
    % featureTable = table(length_skel, thickness, thickness_std, std_ratio, ...
    %                      tortuosity, branch_count, euler_vec, ...
    %                      solidity_vec, aspect_ratio, ...
    %                      'VariableNames', ...
    %                      {'SkelLength','Thickness','ThicknessStd','StdRatio', ...
    %                       'Tortuosity','BranchCount','EulerNum', ...
    %                       'Solidity','AspectRatio'});
    % Only take variables we need for now:
        featureTable = table(thickness, euler_vec, ...
                         'VariableNames', ...
                         {'Thickness','EulerNum'});
    
    %% ════════════════════════════════════════════════════════════════════════
    %  STEP 3  –  Classification
    %
    %  Two modes:
    %   (A) External  – params.AJCatCategories provided  → evaluate in order
    %   (B) Built-in  – fallback using scalar param thresholds
    %
    %  In BOTH modes the last category is an implicit catch-all (else-branch).
    % ════════════════════════════════════════════════════════════════════════
    catAssign = zeros(numObj, 1);   % 0 = unassigned
    
    % if ~isempty(params.AJCatCategories)
    % ── external category definitions ──────────────────────────────────
    numCats = numel(params.AJCatCategories);

    for c = 1 : numCats %- 1        % last cat = catch-all, skip for now
        condCell  = params.AJCatCategories{c};
        catMask   = evalCategoryConditions(featureTable, condCell);

        % Assign only fragments not yet classified
        newAssign = catMask & (catAssign == 0);
        catAssign(newAssign) = c;
    end

    % Last category catches everything remaining
    catAssign(catAssign == 0) = numCats;
    
    % numCats = max(catAssign);
    
    %% ── Label map ────────────────────────────────────────────────────────────
    labelMap = zeros(size(BW));   % double by default
    for i = 1 : numObj
        mask = false(size(BW));
        mask(CC_holey.PixelIdxList{i}) = true;
        labelMap(mask) = catAssign(i);
    end

    %% ── Skeleton-based category fractions ───────────────────────────────────
    % Count skeleton pixels per category; normalise by total skeleton pixels.
    % Using the skeleton -> thick blobs do not outweigh thin strands.
    skelTotal = sum(Skel2(:));
    if skelTotal > 0
        for c = 1 : numCats
            % Skeleton pixels whose label matches this category
            summaryFrac(row, c) = sum(Skel2(labelMap == c)) / skelTotal; % percentages are calculated on number of skeleton pixels that fall in each category
        end
    else
        % No skeleton pixels (tiny/round fragment) – split evenly or leave nan
        summaryFrac(row, :) = nan;
    end
    
    %% ── Visualization ────────────────────────────────────────────────────────
    % % Default colour palette: Cat1=green, Cat2=red, Cat3=blue, then jet for more
    % palette = [0 1 0; 1 0 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 .5 0; .5 0 1];
    % if numCats > size(palette,1)
    %     palette = [palette; jet(numCats - size(palette,1))];
    % end
    % 
    % RGB = zeros([size(BW), 3]);
    % for c = 1 : numCats
    %     mask2d = labelMap == c;
    %     RGB(:,:,1) = RGB(:,:,1) + palette(c,1) * double(mask2d);
    %     RGB(:,:,2) = RGB(:,:,2) + palette(c,2) * double(mask2d);
    %     RGB(:,:,3) = RGB(:,:,3) + palette(c,3) * double(mask2d);
    % end
    
    %% ── Output table ─────────────────────────────────────────────────────────
    % catNames = arrayfun(@(c) sprintf('Cat%d', c), catAssign, 'UniformOutput', false);
    % statsTable = [featureTable, table(catAssign, catNames, ...
    %               'VariableNames', {'CatNum','Class'})];

    maskClasses{row} = labelMap;
end % for row loop

T.MaskClass = maskClasses;
T.FA_summaryFrac = summaryFrac;
fprintf('Done.\n');

end  % ── main function ────────────────────────────────────────────────────


%% ════════════════════════════════════════════════════════════════════════
%  LOCAL: evalCategoryConditions
%
%  Parses one category cell and returns a logical vector (numObj × 1).
%
%  condCell format: flat list of triplets  {feat, op, val, feat, op, val …}
%  All triplets are combined with AND.
%
%  Example:
%    { 'Thickness','<',8, 'StdRatio','<',0.35 }
%    → Thickness < 8  AND  StdRatio < 0.35
% ════════════════════════════════════════════════════════════════════════
function mask = evalCategoryConditions(tbl, condCell)
    n    = height(tbl);
    mask = true(n, 1);

    if mod(numel(condCell), 3) ~= 0
        error('Calc_classifyAdhesions: each category cell must contain triplets {feature, operator, value}.');
    end

    for k = 1 : 3 : numel(condCell)
        feat = condCell{k};
        op   = condCell{k+1};
        val  = condCell{k+2};

        if ~ismember(feat, tbl.Properties.VariableNames)
            error('Calc_classifyAdhesions: unknown feature "%s". Check spelling against statsTable column names.', feat);
        end

        col = tbl.(feat);

        switch op
            case '<',  mask = mask & (col <  val);
            case '<=', mask = mask & (col <= val);
            case '>',  mask = mask & (col >  val);
            case '>=', mask = mask & (col >= val);
            case '==', mask = mask & (col == val);
            case '~=', mask = mask & (col ~= val);
            otherwise
                error('Calc_classifyAdhesions: unsupported operator "%s". Use <  <=  >  >=  ==  ~=', op);
        end
    end
end