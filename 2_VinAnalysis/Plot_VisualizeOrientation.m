function Plot_VisualizeOrientation(T, rowIdx)
    row   = T(rowIdx, :);
    mask  = row.Mask{1};
    angs  = row.Orientations{1};
    muAng = row.MeanAngle;
    AI    = row.AlignIndex;

    figure('Position', [100 100 1100 420]);

    % --- Panel 1: mask ---
    subplot(1,3,1)
    imshow(mask)
    title(sprintf('FOV %d — mask', rowIdx))

    % --- Panel 2: orientation colour map ---
    CC = bwconncomp(mask, 8);
    rp = regionprops(CC, 'Orientation', 'MajorAxisLength', ...
                         'MinorAxisLength', 'Centroid', 'PixelIdxList');
    aspectRatio = [rp.MajorAxisLength] ./ max([rp.MinorAxisLength], 1);
    keep = find(aspectRatio > 1.5);

    rgbMask = repmat(uint8(mask) * 40, 1, 1, 3);   % dark grey background
    cmap    = hsv(180);                              % hue = angle

    subplot(1,3,2); imshow(rgbMask); hold on
    for k = keep
        ang    = mod(rp(k).Orientation, 180);       % 0–179
        colour = cmap(max(1, round(ang)), :);
        % draw centroid dot coloured by orientation
        plot(rp(k).Centroid(1), rp(k).Centroid(2), ...
             'o', 'Color', colour, 'MarkerFaceColor', colour, 'MarkerSize', 5)
    end
    title(sprintf('Orientation map  (AI = %.2f)', AI))
    colormap(gca, hsv); clim([0 180])
    c = colorbar; c.Label.String = 'Angle (°)'; hold off

    % --- Panel 3: rose / polar histogram ---
    subplot(1,3,3)
    polarhistogram(deg2rad(angs * 2), 36, ...   % *2 for axial symmetry
        'Normalization', 'probability', 'FaceColor', [0.2 0.5 0.9])
    title(sprintf('Rose plot  μ = %.0f°  AI = %.2f', muAng, AI))
    % Add mean direction line
    hold on
    maxR = 0.15;
    polarplot([deg2rad(muAng*2) deg2rad(muAng*2)], [0 maxR], ...
              'r-', 'LineWidth', 2)
    polarplot([deg2rad(muAng*2+pi) deg2rad(muAng*2+pi)], [0 maxR], ...
              'r-', 'LineWidth', 2)
    hold off
end