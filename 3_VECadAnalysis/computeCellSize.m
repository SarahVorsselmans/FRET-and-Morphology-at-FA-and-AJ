function CellMask = computeCellSize(B)
% COMPUTECELLMASK
% Returns a binary mask of valid cells (1 = cell interior, 0 = background)
% Keeps only fully enclosed, interior cells (no border-touching cells)

    % --- Defaults
    % if nargin < 2 || isempty(rClose),  rClose  = 10;  end
    % if nargin < 3 || isempty(rThick),  rThick  = 1;  end
    % % if nargin < 4 || isempty(margin),  margin  = 10; end
    % if nargin < 4, opts = struct; end

    rClose = 10;
    rThick = 1;
    borderGapThresh = 40;

    opts.minAreaPx = 1000;
    opts.maxAreaPx = 40000;
    opts.coreInset = 2;
    opts.imgBorderWidth = 2;

    % --- Ensure logical
    if isa(B,'dip_image'), B = dip_array(B); end
    B = B ~= 0;

    [H,W] = size(B);

    % ============================================================
    % --- Step 1: stitch borders at image edges
    % ============================================================
    mask = B;

    mask(1,:)   = stitch1D(mask(1,:),   borderGapThresh);
    mask(end,:) = stitch1D(mask(end,:), borderGapThresh);
    mask(:,1)   = stitch1D(mask(:,1),   borderGapThresh);
    mask(:,end) = stitch1D(mask(:,end), borderGapThresh);

    % Regularize borders
    Bc = imclose(mask, strel('disk', rClose));
    Bc = bwmorph(Bc,'bridge');
    Bt = imdilate(Bc, strel('disk', rThick));

    % ============================================================
    % --- Step 2: find enclosed regions (cells)
    % ============================================================
    holes = imfill(Bt,'holes') & ~Bt;
    holes = bwareaopen(holes, opts.minAreaPx);

    % ============================================================
    % --- Step 3: define valid interior region (core)
    % ============================================================
    % core = true(H,W);
    % 
    % margin = min([margin, floor(H/4), floor(W/4)]);
    % if margin > 0
    %     core(1:margin,:) = false;
    %     core(end-margin+1:end,:) = false;
    %     core(:,1:margin) = false;
    %     core(:,end-margin+1:end) = false;
    % end

    % if opts.coreInset > 0
    %     innerCore = imerode(core, strel('disk', opts.coreInset));
    % else
    %     innerCore = core;
    % end

    % Image border mask
    % imgBorder = false(H,W);
    % w = max(1, round(opts.imgBorderWidth));
    % imgBorder(1:w,:) = true;
    % imgBorder(end-w+1:end,:) = true;
    % imgBorder(:,1:w) = true;
    % imgBorder(:,end-w+1:end) = true;

    % ============================================================
    % --- Step 4: filter valid cells (compact version)
    % ============================================================
    CC = bwconncomp(holes);

    keepMask = false(H,W);

    for k = 1:CC.NumObjects
        idx = CC.PixelIdxList{k};
        Ak  = numel(idx);

        if Ak < opts.minAreaPx || Ak > opts.maxAreaPx
            continue;
        end

        % if any(~innerCore(idx))
        %     continue;
        % end

        % if any(imgBorder(idx))
        %     continue;
        % end

        keepMask(idx) = true;
    end

    CellMask = keepMask;
end


% ============================================================
% --- stitch1D helper (unchanged)
% ============================================================
function v = stitch1D(v, maxGap)
    v = v(:)';
    n = numel(v);
    idx = 1;

    while idx <= n
        if v(idx) == 0
            runStart = idx;
            while idx <= n && v(idx) == 0
                idx = idx + 1;
            end
            runEnd = idx - 1;
            runLen = runEnd - runStart + 1;

            leftFlanked  = runStart > 1 && v(runStart - 1) == 1;
            rightFlanked = runEnd   < n && v(runEnd + 1) == 1;

            if leftFlanked && rightFlanked && runLen <= maxGap
                v(runStart:runEnd) = 1;
            end
        else
            idx = idx + 1;
        end
    end
end


% ============================================================
% --- defaults helper
% ============================================================
function s = setDefaults(s, def)
    f = fieldnames(def);
    for k = 1:numel(f)
        if ~isfield(s,f{k}) || isempty(s.(f{k}))
            s.(f{k}) = def.(f{k});
        end
    end
end

% function out = computeCellSize(B, px, rClose, rThick, margin, I0, showFig, opts)
% % QUICKAVERAGEFROMCLOSEDBORDERS
% % Average cell area from VE-cadherin borders using "enclosed holes".
% % Counts only cell interiors that are fully enclosed by borders and fully
% % inside a core region. Adds a border-support check to exclude weak/ambiguous regions.
% %
% % Usage:
% %   out = quickAverageFromClosedBorders(B, px)
% %   out = quickAverageFromClosedBorders(B, px, rClose, rThick, margin, I0, showFig, opts)
% %
% % Inputs
% %   B        : border mask (true where VE-cadherin borders are). logical or dip_image.
% %   px       : µm/pixel (optional; [] -> results in px²)
% %   rClose   : gap closing radius (default 2)
% %   rThick   : border thickening radius (default 1)
% %   margin   : outer margin for core (default 10 px)
% %   I0       : (optional) grayscale VE-cadherin image for the overlay background
% %   showFig  : (optional) true/false to display figure (default true)
% %   opts     : struct with optional fields:
% %                .minAreaPx     (default 100)   - min cell area
% %                .maxAreaPx     (default inf)   - max cell area
% %                .coreInset     (default 2)     - erode core to ensure full inclusion
% %                .dEnclose      (default 2)     - distance (px) for border support
% %                .fracEnclose   (default 0.80)  - min fraction of perimeter close to borders
% %                .bridge        (default true)  - apply bwmorph(...,'bridge')
% %
% % Output struct 'out'
% %   mean_px, median_px, areas_px, (and *_um2 if px given)
% %   T                 : per-cell table
% %   L                 : final label image (kept cells)
% %   core, innerCore   : masks
% %   overlayRGB        : RGB overlay image
% %   fig               : figure handle (if shown)
% %   counts            : struct with kept/dropped statistics
% 
%     % --- Defaults
%     if nargin < 2, px = []; end
%     if nargin < 3 || isempty(rClose),  rClose  = 2;  end
%     if nargin < 4 || isempty(rThick),  rThick  = 1;  end
%     if nargin < 5 || isempty(margin),  margin  = 10; end
%     if nargin < 6, I0 = []; end
%     if nargin < 7 || isempty(showFig), showFig = true; end
%     if nargin < 8, opts = struct; end
%     % opts = setDefaults(opts, struct( ...
%     %     'minAreaPx',100, 'maxAreaPx',inf, 'coreInset',2, ...
%     %     'dEnclose',2, 'fracEnclose',0.80, 'bridge',true));
%     opts = setDefaults(opts, struct( ...
%     'minAreaPx',100, 'maxAreaPx',inf, 'coreInset',2, ...
%     'dEnclose',2, 'fracEnclose',0.70, 'bridge',true, ...
%     'minCoreCoverage',0.90, ...      % NEW: allow small core crossing
%     'maxBorderTouchFrac',0.07, ...   % NEW: allow slight image-border touch
%     'imgBorderWidth',2,...           % NEW: how thick to consider the image border
%     'borderGapThresh',40));          % Gap at FOV border to be filled
% 
%     % --- DIPimage compatibility
%     if isa(B,'dip_image'), B = dip_array(B); end
%     B = B ~= 0;
%     if ~isempty(I0) && isa(I0,'dip_image'), I0 = dip_array(I0); end
% 
%     % --- Step 0: close FOV border gaps 
%     B(1,:)   = stitch1D(B(1,:),   opts.borderGapThresh);   % top
%     B(end,:) = stitch1D(B(end,:), opts.borderGapThresh);   % bottom
%     B(:,1)   = stitch1D(B(:,1),   opts.borderGapThresh);   % left
%     B(:,end) = stitch1D(B(:,end), opts.borderGapThresh);   % right
% 
%     % --- Step 1: strengthen/regularize border network
%     Bc = imclose(B, strel('disk', rClose));
%     if opts.bridge
%         Bc = bwmorph(Bc,'bridge');  % connect diagonal near-misses
%     end
%     Bt = imdilate(Bc, strel('disk', rThick)); % thicken a bit to seal thin gaps
% 
%     % --- Step 2: candidate cells = "holes" inside thickened border network
%     % Holes are areas completely enclosed by borders. Open areas aren't holes.
%     holes = imfill(Bt,'holes') & ~Bt;
%     holes = bwareaopen(holes, opts.minAreaPx); % drop tiny specks
% 
%     % --- Step 3: core and inner-core (to keep only fully contained cells)
%     [H,W] = size(B);
%     imgBorder = false(H,W);
%     w = max(1, round(opts.imgBorderWidth));
%     imgBorder(1:w,:) = true; imgBorder(end-w+1:end,:) = true;
%     imgBorder(:,1:w) = true; imgBorder(:,end-w+1:end) = true;
%     %[H,W] = size(B);
%     core = true(H,W);
%     margin = min([margin, floor(H/4), floor(W/4)]);
%     if margin > 0
%         core(1:margin,:) = false; core(end-margin+1:end,:) = false;
%         core(:,1:margin) = false; core(:,end-margin+1:end) = false;
%     end
%     if opts.coreInset > 0
%         innerCore = imerode(core, strel('disk', opts.coreInset));
%     else
%         innerCore = core;
%     end
% 
%     % --- Step 4: perimeter border-support filter
%     Lc = bwconncomp(holes);
%     statsAll = regionprops(Lc, 'Area','PixelIdxList');
% 
%     keep   = false(numel(statsAll),1);
%     reason = strings(numel(statsAll),1);
% 
%     for k = 1:numel(statsAll)
%         Pk = statsAll(k).PixelIdxList;
%         Ak = statsAll(k).Area;
% 
%         % --- Area filter
%         if Ak < opts.minAreaPx || Ak > opts.maxAreaPx
%             reason(k) = "area";
%             continue;
%         end
% 
%         % --- Must be fully inside inner core
%         if any(~innerCore(Pk))
%             reason(k) = "core";
%             continue;
%         end
% 
%         % --- Reject anything touching image border
%         if any(imgBorder(Pk))
%             reason(k) = "border";
%             continue;
%         end
% 
%         % --- Passed all filters
%         keep(k) = true;
%         reason(k) = "kept";
%     end
% 
%     keptLabels = find(keep);
%     droppedLabels = find(~keep);
%     maskKept = ismember(Lc, keptLabels);
%     L = bwlabel(maskKept);
% 
%     % --- Measurements on kept cells
%     stats = regionprops(L, 'Area','Centroid','BoundingBox');
%     areas_px = [stats.Area]';
%     out.mean_px = mean(areas_px);
%     out.median_px = median(areas_px);
%     out.areas_px = areas_px;
%     if ~isempty(px)
%         out.mean_um2 = out.mean_px * px^2;
%         out.median_um2 = out.median_px * px^2;
%         out.areas_um2 = areas_px * px^2;
%     end
% 
%     % Per-cell table
%     Label = (1:numel(stats))';
%     Area_px = areas_px;
%     Centroid = reshape([stats.Centroid],2,[])';
%     BBox = vertcat(stats.BoundingBox);
%     T = table(Label, Area_px, Centroid, BBox);
%     if ~isempty(px), T.Area_um2 = Area_px * px^2; end
%     out.T = T;
% 
%     % --- Overlays
%     base = prepareBaseForOverlay(I0, B);     % grayscale base % holes
%     overlayRGB = makeOverlay(base, L, innerCore, 0.28);
%     out.overlayRGB = overlayRGB;
%     %out.L = L; out.core = core; out.innerCore = innerCore;
% 
%     % --- Counts summary
%     out.counts = struct( ...
%     'candidates', numel(statsAll), ...
%     'kept', numel(keptLabels), ...
%     'dropped_total', numel(droppedLabels));
% 
%     if showFig
%         out.fig = figure('Color','w','Name','Cells: enclosed-hole selection','Units','normalized','Position',[0.07 0.09 0.86 0.74]);
%         tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
% 
%         % Left: borders & cores
%         ax1 = nexttile; imshow(base,'Parent',ax1); hold(ax1,'on');
%         visBoundaries(ax1, Bt, [0 0 0]); % thickened borders in black
%         rectangle(ax1,'Position',[1+margin,1+margin,W-2*margin,H-2*margin], ...
%                   'EdgeColor',[0 1 1],'LineWidth',1.5,'LineStyle','--');
%         rectangle(ax1,'Position',[1+margin+opts.coreInset,1+margin+opts.coreInset, ...
%                   W-2*(margin+opts.coreInset), H-2*(margin+opts.coreInset)], ...
%                   'EdgeColor',[0 0.6 1],'LineWidth',1.0,'LineStyle',':');
%         title(ax1, sprintf('Borders + Core (margin=%d, inset=%d)', margin, opts.coreInset), 'Interpreter','none');
% 
%         % Right: kept vs dropped
%         ax2 = nexttile; imshow(overlayRGB,'Parent',ax2); hold(ax2,'on');
%         visBoundaries(ax2, L>0, [1 0 0]); % kept in red
%         % Outline dropped candidates in gray
%         if ~isempty(droppedLabels)
%             Ldrop = ismember(Lc, droppedLabels);
%             visBoundaries(ax2, Ldrop, [0.6 0.6 0.6]);
%         end
%         % Labels with areas
%         for k = 1:numel(stats)
%             c = stats(k).Centroid;
%             if isempty(px)
%                 labelStr = sprintf('%d (%.0f px^2)', k, stats(k).Area);
%             else
%                 labelStr = sprintf('%d (%.0f µm^2)', k, stats(k).Area*px^2);
%             end
%             text(ax2, c(1), c(2), labelStr, 'Color','w', 'FontSize',10, 'FontWeight','bold', ...
%                 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
%                 'BackgroundColor',[0 0 0 0.35]);
%         end
%         title(ax2, sprintf('Kept N = %d  |  Dropped = %d', numel(stats), numel(droppedLabels)));
%         annotation(out.fig,'textbox',[0.32 0.91 0.36 0.07],'String', ...
%             'Red = kept cells; Gray outline = dropped candidates; Cyan dashed = core; Cyan dotted = inner core', ...
%             'FitBoxToText','on','EdgeColor','none','HorizontalAlignment','center');
% 
%         % === Diagnostic labels for ALL candidates (kept + dropped) ===
%         if showFig
%             % After you've created ax2 and displayed the right image:
%             statsAll2 = regionprops(Lc, 'Centroid');
%             Db = bwdist(Bt); % already computed above
% 
%             for k = 1:numel(statsAll)
%                 c = statsAll2(k).Centroid;
% 
%                 % compute fracEnclose again for display (or store it during the loop)
%                 Rk = false(size(B)); Rk(statsAll(k).PixelIdxList) = true;
%                 perim = bwperim(Rk); if ~any(perim(:)), perim = Rk; end
%                 frac = mean(Db(perim) <= opts.dEnclose);
% 
%                 col = [1 1 0]; % yellow text for diagnostics
%                 str = sprintf('#%d  %s  (f=%.2f)', k, reason(k), frac);
%                 text(ax2, c(1), c(2)-12, str, 'Color', col, 'FontSize',9, ...
%                      'FontWeight','bold', 'HorizontalAlignment','center', ...
%                      'BackgroundColor',[0 0 0 0.25], 'Margin',1);
%             end
%         end
%     end
% end
% 
% % ---------- Helpers (same as before) ----------
% 
% function base = prepareBaseForOverlay(I0, Ifallback)
%     if ~isempty(I0)
%         I0 = im2double(I0);
%         if ndims(I0) == 3, I0 = rgb2gray(I0); end
%         base = mat2gray(I0);
%     else
%         base = im2double(Ifallback); % falls back to holes or interior mask
%     end
% end
% 
% function RGB = makeOverlay(baseGray, L, core, alpha)
%     baseGray = mat2gray(baseGray);
%     baseRGB  = repmat(baseGray, [1 1 3]);
%     C = label2rgb(L,'lines','k','shuffle'); C = im2double(C);
%     mask = (L > 0) & core;
%     RGB = baseRGB;
%     for ch=1:3
%         chB = RGB(:,:,ch); chC = C(:,:,ch);
%         chB(mask) = (1 - alpha)*chB(mask) + alpha*chC(mask);
%         RGB(:,:,ch) = chB;
%     end
% end
% 
% function visBoundaries(ax, BW, colorRGB)
%     P = bwperim(BW);
%     [y,x] = find(P);
%     plot(ax, x, y, '.', 'Color', colorRGB, 'MarkerSize', 1.5);
% end
% 
% function s = setDefaults(s, def)
%     f = fieldnames(def);
%     for k = 1:numel(f)
%         if ~isfield(s,f{k}) || isempty(s.(f{k})), s.(f{k}) = def.(f{k}); end
%     end
% end
% 
% 
% function v = stitch1D(v, maxGap)
% % STITCH1D  Fill runs of zeros in a 1-D logical vector that are flanked on
% %           both sides by ones, provided the run length <= maxGap.
% %
% %   v      - 1-D logical vector (row or column, both handled)
% %   maxGap - maximum gap length to close
% 
% v    = v(:)';          % work as row vector
% n    = numel(v);
% idx  = 1;
%     while idx <= n
%         if v(idx) == 0
%             % Find end of this zero-run
%             runStart = idx;
%             while idx <= n && v(idx) == 0
%                 idx = idx + 1;
%             end
%             runEnd = idx - 1;
%             runLen = runEnd - runStart + 1;
% 
%             % Stitch only if flanked on both sides AND within threshold
%             leftFlanked  = runStart > 1   && v(runStart - 1) == 1;
%             rightFlanked = runEnd   < n   && v(runEnd   + 1) == 1;
% 
%             if leftFlanked && rightFlanked && runLen <= maxGap
%                 v(runStart:runEnd) = 1;
%             end
%         else
%             idx = idx + 1;
%         end
%     end
% end