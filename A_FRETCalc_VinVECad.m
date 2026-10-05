function [Int, FRET_Index_Data, mask, outputMetrics, FA] = A_FRETCalc_VinVECad(data, flags)

disp(['--- Analyzing ' data.sampleIndex ' ---'])
%To collapse all the sections press these three keys; Ctrl+Shift+N hehe 
%% 1. The Beginning
%addpath("functions\");
load('BlackJet.mat');
color.blueTransp = [0.7 0.7 0.9 0.5];
color.orangeTransp = [0.9 0.7 0.7 0.5];
IntTag = 'ch1_ch0';
GTag = 'ch1_ch1';
STag = 'ch1_ch2';

%% User Inputs
Intensity_Threshold = 8; %just to remove very dim pixels % was 8
user.dataPath = data.dataPath;
fileTag = '.tif';
user.outputPath = data.outputPath;
user.outputPath = SetOut(user.outputPath, user.dataPath);
user.outputType = '.tiff';
user.FRET_Threshold = 50;
user.FRET_ColorbarLim = 40;
user.Image_Resolution = 512;
user.Image_Size = 123.26; 
% user.pixelsize = 0.361; % in µm 
user.medfiltValue = 9; 
%user.minimumFA = 25;  % %CHANGE FOR VE: 25, Vin: 8 (now 15)
user.minSize_Vin = 8;
user.minSize_VE = 25;

%% Size
if contains(data.ConditionFolder, 'Vin')
    user.minimumFA = user.minSize_Vin;
else
    user.minimumFA = user.minSize_VE;
end

%% Extracting Data
%addpath("functions\FocalAdhesionMasking");

% for Loop_FDA

Int = data.int;
FLIMG = data.FLIMG;
FLIMS = data.FLIMS;

user.saveName = data.name; %(to save the name of the sample)

if data.flags.AcceptorData.switch
    acceptor_channel = data.acceptor;
else
end

%% 2. Acceptor channel analysis

if data.flags.AcceptorData.switch
    %addpath("functions\AcceptorAnalysis");
    acceptor_channel = ShiftCorrChannel(Int, acceptor_channel);
else
end

% Define which data is used to ROI on
if flags.ROIDorA == 0
    Input = acceptor_channel;
else
    Input = Int;
end

% for the moment: check that shift correction is working properly
% by uncommenting the 'joinchannels' functions in 'ShiftCorrChannel'

%% 3. Remove/only keep regions
% if ismember(flags.manualROI, [1, 2]) %- version1
%     RegionMask = interactive_roi(Int); % Opens function to manually draw ROI. The ROI will be 1.
%     if flags.manualROI == 2
%         RegionMask = ~RegionMask; % Inverts mask, so that ROI is 0
%     else
%     end
% Assume 'flags' is defined and ROIFile is a string with the file path

% If bg = black
%data.mask = ~data.mask; %do not activate

[RegionMask] = ROIandSeg(data, flags, user.minimumFA);

%% 4. Cell Segment

if flags.CellSegment == 1
    if data.CellSegThreshSwitch == 1
        [row,~] = find(data.CellSegThresh == data.index);
        if ~isempty(row)
            threshFactor = data.CellSegThresh(row, 2);
            disp(['• Using CellSegThresh factor from .txt file: ' num2str(threshFactor)])
            CellMask = CellSegment(Int, threshFactor);
        else
            threshFactor = NaN;
            %CellSegThesh is there, but either this cell is missing or the threshold value is missing. 
            disp('• This cell is missing in CellSegThresh file')
            if flags.AutoCellSegment == 1
                CellMask = CellSegment(Int, threshFactor);
            elseif flags.AutoCellSegment == 2
                path = [user.dataPath '\' data.sampleIndex '_CellSegment.tiff'];
                CellMask = imread(path);
                CellMask =~ CellMask;
            else
            end
        end
    else
        %if there is no CellSegThresh file at all. 
        disp('• No CellSegThresh.txt found, using AutoCellSegment')
        if flags.AutoCellSegment == 1
            CellMask = CellSegment(Int, NaN);
        elseif flags.AutoCellSegment == 2
            path = [user.dataPath '\' data.sampleIndex '_CellSegment.tiff'];
            CellMask = imread(path);
            CellMask = ~CellMask;
        else
            error('No segmentation method available (no txt file, AutoCellSegment off)')
        end
    end
else
    CellMask = true(512,512); % Make a CellMask that masks everything (nothing is excluded)
end

%% Creating the FA Mask
if flags.AutoFASegments == 1 % automatic segmentation - still artificially closes ROIs so be careful
    
    [bw_FA_mask1] = FASegment(Input, CellMask, 1);%takes the brighest
    [bw_FA_mask2] = FASegment(Input, CellMask, 0.5);%takes the lowest
    [bw_FA_mask3] = FASegment(Input, CellMask, 1.2);%takes in between

    mask_rim = logical(bw_FA_mask2 - bw_FA_mask1);

    skel = bwskel(mask_rim);
   
    skel(bw_FA_mask1 == 1) = 0;

    FAmask1 = bw_FA_mask3;
    FAmask1(skel == 1) = 0;
 
    CC = bwconncomp(FAmask1, 4);
    L = labelmatrix(CC);
    stats = regionprops(CC,'area');
    BW = ismember(L, find([stats.Area] > 0  ));  %& [stats.Area] < medArea*60)
    CC = bwconncomp(BW, 4);
    L = labelmatrix(CC);
    
    outline = logical(imdilate(BW, strel('disk', 1)) - BW);
 
    CCoutline = bwconncomp(outline, 8);
    Lout = labelmatrix(CCoutline);
    stats = regionprops(CCoutline,'area');
  
    BW2 = BW;
    BW2(bw_FA_mask1 == 1) = 1;
    BW2(outline == 1) = 0;

    % close holes enclosed in FA area
    se = strel('disk', 1); % Structural element for dilation/erosion
    BW3 = imdilate(BW2, se); % Opening operation
    BW3 = imerode(BW3, se); % Opening operation
    BW3(outline == 1) = 0;
    CC = bwconncomp(BW3, 4);
    L3 = labelmatrix(CC);

    % remove FAs that are too small
    stats = regionprops(CC,'area');
    BW3 = ismember(L3, find([stats.Area] > user.minimumFA  ));  %& [stats.Area] < medArea*60)

    CC2 = bwconncomp(BW3, 4);
    L2 = labelmatrix(CC2);
    %   dipshow(L2); %use this if you want to adjust minimum FA threshold :)
    %finalValuestring = 'Adaptive';
    finalValuestring = num2str(Intensity_Threshold);
    % the idea: to remove skeleton when not needed
    empty = zeros(CC2.ImageSize);
    summedBW = empty;
    for i = 1:CC2.NumObjects
        singleBW = empty;
        singleBW(cell2mat(CC2.PixelIdxList(i))) = 1;
        se = strel('disk', 1); % Structural element for dilation/erosion
        singleBW = imdilate(singleBW, se); % Opening operation
        summedBW = summedBW + singleBW; % overlapping areas will have a value greater than 1
    end
    skel2 = empty;
    skel2(summedBW > 1) = 1;
    newBW = bw_FA_mask3 - skel2;
    newBW(outline == 1) = 0;
    % filter based on size
    CCnew = bwconncomp(newBW, 4);
    Lnew = labelmatrix(CCnew);
    stats = regionprops(CCnew,'area');
    newBW = ismember(Lnew, find([stats.Area] > user.minimumFA )); % remove small FAs
    % smallBW = ismember(Lnew, find([stats.Area] == 1 )); % remove small FAs
    % newBW(smallBW == 1) = 0;
    CCnew2 = bwconncomp(newBW, 4);
    Lnew2 = labelmatrix(CCnew2);
    dipshow(Lnew2)
    % dipshow(skel2);
    FASegmentoutput = newBW; % output this to the rest of the code

elseif flags.AutoFASegments == 0 % use  slider

    % A Slider user interface to find a desirable Intensity Threshold
    finalValue = SliderImagePlot2(Input);
    Intensity_Threshold = finalValue;

    finalValuestring = num2str(finalValue); %just to write the threshold in the name of the FRET Index Data file

elseif flags.AutoFASegments == 2 % use Ilastik generated .tiff files % BUT IS DOESN'T USE ANY ILASTIK FILE???

    CC = bwconncomp(CellMask); %make a list of all the objects in your mask
    stats = regionprops(CC,'area');%gives you the area of each objective

    [areas] = deal([stats.Area]); %we dont know what does 'deal' do
    [area,indx] = max(areas); %just tells you which object has the maximum area, and you get the index information of this object

    %bw2 = false(size(CellMask)); %you create an empty mask (everything is false or Zero)
    %bw2(CC.PixelIdxList{indx}) = true; %then you make only the object with maximum area 1 ('true')
    %SE = strel('disk',0);%creates a shape, but now its set to Zero
    %bw2 = imdilate(bw2,SE);%applies the shape everywhere where the index is 1 % this does nothing by itself?
    %(so your largest objective basically)
    %CellMask = imfill(bw2,'holes');%any closed shape with hole will be filled
    %CellMask = bw2;

    finalValuestring = num2str(Intensity_Threshold);

else
end

%% 6. Applying the Mask 

if flags.AutoFASegments == 1 %automatic segmentation
    FAmask = FASegmentoutput;
    if ismember(flags.manualROI, [1, 2]) % NEW
        FAmask = FAmask & ~RegionMask; %(| for both - & for overlapping - , for non-verlapping - &~ for subtraction)
    else
    end
    mask = FAmask;
else
    extMask = flags.useExternalMask;
    if  extMask == false
        mask = Int;
        mask(mask < 1) = 0;
        mask(mask <= Intensity_Threshold) = 0;
        mask(mask >= Intensity_Threshold) = 1;
        mask = bwareaopen(mask, user.minimumFA);

    elseif flags.AutoFASegments == 0 %slider function

        mask = data.mask;
        mask = mask(:,:,1);
        mask(mask == 255) = 1;
        mask = cast(mask, 'logical');
        mask = ~mask;
        mask = imfill(mask,"holes");
        mask(Int <= Intensity_Threshold) = 0; % still applies the threshold after the Ilastik mask
        mask = bwareaopen(mask, user.minimumFA);
        disp(['Removing all FAs <' num2str(user.minimumFA) ' [pixels²] from mask'])

        fused = imfuse(mask, Int);
        figure("Name", 'Mask on Intensity');
        imshow(fused)
        axis square
        hold on

    elseif flags.AutoFASegments == 2 %ilastik mask
        mask = data.mask;
        mask = mask(:,:,1);
        mask(mask == 255) = 1;
        mask = cast(mask, 'logical'); % converts the mask to logical (true/false) instead of numeric 0/1.
        %mask = ~mask; % inverts the mask
        %mask = imfill(mask,"holes");
        mask(Int <= Intensity_Threshold) = 0;
        mask = bwareaopen(mask, user.minimumFA, 8); % removes small objects % 8 = also diagonal % already done in ROIandSeg

        % For Vin: get rid of outer pixels of a FA:
        % BW_outline = bwperim(mask, 4);  % one-pixel outline; diagonal-only contact is ignored
        % mask = mask &~ BW_outline;
    end

    FAmask = mask; %just changing the name because the following lines uses INTmask instead
    if ismember(flags.manualROI, [1, 2]) % NEW
        FAmask = FAmask & RegionMask; %(| for both - & for overlapping - , for non-verlapping - &~ for subtraction)
    else
    end
end
faList = bwconncomp(FAmask);%only after the mask is applied, we extract the list of FAs :)

if flags.flagSaveMasks == 1
    fig000 = figure;
    imagesc(FAmask);
    axis square
    colormap(gray);
    axis off;
    hold on;
    saveas(fig000, [user.outputPath filesep user.saveName '_1. FAMask' user.outputType])
else
end

%% 7. Cell mask image
if flags.CellSegment == 1
    CellPixArea = sum(CellMask(:) == 1);
    
    fig00 = figure;
    imagesc(CellMask);
    axis square
    colormap(gray);
    axis off;
    hold on;

    if flags.flagSaveMasks == 1
    saveas(fig00, [user.outputPath filesep user.saveName '_0. Cell Mask' user.outputType])
    else
    end
else
end

%% 8. Cytosol mask

CytosolMask = CellMask;
CytosolMask(mask == 1) = 0;
CytosolMaskbefore = CytosolMask; % needed for cytosol intensity (including low background signal)

%% Cytosol remove low intensity
img = rescale(Int);
img = medfilt2(img);
img(CellMask == 0) = NaN;
cytosolThreshold = graythresh(img);
cytosolHighMask = imbinarize(Int, cytosolThreshold*0.0025);
cytosolHighMask = ~cytosolHighMask;
cytosolHighMask(FAmask == 1) = 0;
cytosolHighMask(CellMask == 0) = 0;
dipshow(cytosolHighMask)


if flags.CytosolClean
    CytosolMask(cytosolHighMask == 1) = 0;
else
end

% cytosol
if flags.flagSaveMasks == 1
    fig0000 = figure;
    imagesc(CytosolMask);
    axis square
    colormap(gray);
    axis off;
    hold on;
    saveas(fig0000, [user.outputPath filesep user.saveName '_2. Cytosol Mask' user.outputType])
else
end

%% 9. Phasor Calculation

% Median filter for the G/S coordinates
FLIMS(FAmask == 0) = NaN;
FLIMG(FAmask == 0) = NaN; 
FLIMS = nanmedfilt2(FLIMS,[9 9]);
FLIMG = nanmedfilt2(FLIMG,[9 9]);
FLIMS(FAmask == 0) = NaN;
FLIMG(FAmask == 0) = NaN; 
% reshape the pixel coordinates into a list & convert to phasor space
size1 = size(FLIMG,1); %simply measures the size of FLIMG (512), which is obviously same for FLIMS.
ListG = (reshape(FLIMG,size1*size1,1)/65535)*2-1; % to convert the raw data to phasor space
ListS = (reshape(FLIMS,size1*size1,1)/65535)*2-1; % to convert the raw data to phasor space

PhasorResolution = 50; %only for visual purposes. 
Laserrep = 19500000; % Laser rep in Hz, 19500000 for 20 MHz, 39000000 for 40 MHz, 78000000 for 80 MHz 
BGg = 0.5; BGs = 0.25;
DonorPureLifetime = 3.0; %IMPORTANT % zero-FRET limit
DonorG =  1/(1+((2*pi*Laserrep)^2)*((DonorPureLifetime*(10^(-9)))^2)); % G = 1/(1+(ωτ)^2) with ω = 2πf
DonorS =  ((2*pi*Laserrep)*(DonorPureLifetime*(10^(-9)))/(1+((2*pi*Laserrep)^2)*((DonorPureLifetime*(10^(-9)))^2))); % S = ωτ/(1+(ωτ)^2) with ω = 2πf
PercentageBG = 1.5; %IMPORTANT, use 1.5 for 20 MHz, 5 for 40 MHz and 7 for 80 MHz. % 0.5??
PercentageDO = 0;

% G- and S-coordinates of FRET trajectory on the semicircle (ideal case) 
TraG=1./(1+(0:0.001:1).^2*(1/DonorG-1));
TraS=(0:0.001:1)*sqrt((1/DonorG)-1)./(1+(0:0.001:1).^2*(1/DonorG-1));

% G- and S-coordinates of the actual FRET trajectory (in the case of BG) 
xT=((100-PercentageBG-PercentageDO)*(0:0.001:1).*TraG +(PercentageBG+PercentageDO)*BGg)./((100-PercentageBG-PercentageDO)*(0:0.001:1)+(PercentageBG+PercentageDO));
yT=((100-PercentageBG-PercentageDO)*(0:0.001:1).*TraS +(PercentageBG+PercentageDO)*BGs)./((100-PercentageBG-PercentageDO)*(0:0.001:1)+(PercentageBG+PercentageDO));

%% FRET Calculation
[T,dist] = dsearchn([xT(1:1001); yT(1:1001)]', [ListG(:,1), ListS(:,1)]); %returns !!indices!! of closest points of trajectory line to every g&s coordinate of the actual sample
Trajectory_Image = reshape(T, size(Int, 1), size(Int, 2)); %what you have is an image of each pixel with an assigned 'cell index' from xT and yT. NOT assigned g and s coordinates.
FRET_Index_Data = abs(Trajectory_Image-1001)/10;

%FRET_Index_Data is the raw data, but remember that everything outside the mask is 100% FRET, so we need to remove them. 

FRET_image = FRET_Index_Data; 

[FRET_image_FA, FRETav, FRET_list, Title_FA, Plow, Phigh] = FRETcalculation(FRET_image, mask, user);

[FRET_image_cyto, FRETav_cyto,  Title_cyto, ~, ~] = FRETcalculation(FRET_image, CytosolMask, user);

[FRET_image_bg, FRETav_bg,  Title_bg, ~, ~] = FRETcalculation(FRET_image, ~CellMask, user);

% % Now we save the FRET_Index_Data (unmasked!!!) to use for TFM-FLIM code. 
% if data.flags.AcceptorData.switch
%    save([user.outputPath filesep user.saveName '_size_' num2str(user.Image_Size) '_threshold_' finalValuestring '.mat'], "FRET_Index_Data", "acceptor_channel", "Int");
% else
%    save([user.outputPath filesep user.saveName '_size_' num2str(user.Image_Size) '_threshold_' finalValuestring '.mat'], "FRET_Index_Data", "Int");
% end


%% Phasor and Histogram Images
if flags.flagPhasorHistogram
    fig = figure;
    histogram2(ListG(:,1),ListS(:,1),PhasorResolution,'DisplayStyle','tile','EdgeColor','none','FaceColor','flat')
        xlabel('G');
        ylabel('S');
        xlim([0 1]);
        ylim([0 1]);
        axis square
        colormap jet;
        hold on
        plot(xT(1:1001),yT(1:1001));
        hold on
        plot(BGg,BGs,'r*');
        hold on
        plot(DonorG,DonorS,'r*');
        hold on
        line([xT(1001) BGg], [yT(1001) BGs]); %take the last element of trajectory, (this is the 0%/beginning of trajectory) and connect with mixture of DO and BG (end of trajectory)
        hold on
        y=0:0.001:1; z=0:0.001:1; % = Semicircle
        plot(y,sqrt(z-z.^2),'k','LineStyle','- -','Linewidth',2, 'color', 'k');
        Xlim = [0 1];
        Ylim = [0 1];
        set(0, 'DefaultFigureRenderer', 'painters');
    hold on;

    saveas(fig, [user.outputPath filesep user.saveName '_Phasor' user.outputType])

    % Extract FRET efficiency values inside FA mask
    FRET_values = FRET_image_FA(~isnan(FRET_image_FA));
    
    % ----- Create histogram and get numeric data -----
    numBins = 100;
    [FRET_counts, FRET_edges] = histcounts(FRET_values, numBins);
    binCenters = FRET_edges(1:end-1) + diff(FRET_edges)/2;
    FRET_counts_norm = FRET_counts / sum(FRET_counts);
    
    
    % ----- Plot histogram -----
    fig = figure;
    bar(binCenters, FRET_counts, 'FaceColor',[0.3 0.3 0.8], 'EdgeColor','none'); % nicer look
    xlabel('FRET (%)');
    ylabel('Pixels (Count)');
    xlim([0 40]);
    axis square;
    
    % outputMetrics.Histogram_BinCenters = binCenters;
    % outputMetrics.Histogram_Counts = FRET_counts;
    % outputMetrics.Histogram_NormalizedCounts = FRET_counts_norm;

    saveas(fig, [user.outputPath filesep user.saveName '_Histogram' user.outputType]);
end


%% 10. Fig 0 Intensity Image

A = Int;
fig = figure;

low = prctile(A(:), 1);
high = prctile(A(:), 99);

imagesc(A, [low high]);
colormap(gray);
colorbar;
axis square;
axis off;

% imagesc(A);
% axis square;
% colorbar;
% colormap (gray);
% axis off;
% hold on;

saveas(fig, [user.outputPath filesep user.saveName '_3. Intensity' user.outputType])


%% Figure 1. FRET/Intensity Overlay (adapted min/max)
if flags.flagFRETimage
    scaled_INT = double(Int)./(double(max(max(Int))));
    A = scaled_INT*2; % 2
    fig6 = figure;
    imagesc(A);
    title(Title_FA);
    axis square
    colormap(gray);
    axis off;
    hold on;
    h1 = imagesc(FRET_image_FA);
    colormap(BlackJet);
    c = colorbar;
    caxis([Plow Phigh]);
    c.Label.String = 'FRET Efficiency %';
    axis square;
    axis off;
    hold on;
    set(h1, 'AlphaData', A);
    hold off
        saveas(fig6, [user.outputPath filesep user.saveName '_4.1 FRET Intensity Overlay Adapted Range' user.outputType])
    else
end
%% Figure 1.2 Just FRET (adapted min/max)

if flags.flagFRETimage
    fig66 = figure;
    h2 = imagesc(FRET_image_FA);
    title(Title_FA);
    axis square
    colormap(BlackJet);
    c = colorbar;
    caxis([Plow Phigh]);
    c.Label.String = 'FRET Efficiency %';
    axis off;
    hold on;
    hold off

        saveas(fig66, [user.outputPath filesep user.saveName '_4.2 FRET Adapted Range' user.outputType])
    else
end

%% Figure 2. FRET/Intensity Overlay (fixed min/max)
if flags.flagFRETimage
    scaled_INT = double(Int)./(double(max(max(Int))));
    A = scaled_INT*2;
    fig5 = figure;
    imagesc(A);
    title(Title_FA);
    axis square
    colormap(gray);
    axis off;
    hold on;
    h = imagesc(FRET_image_FA);
    colormap(BlackJet);
    c = colorbar;
    caxis([0 user.FRET_ColorbarLim]);
    c.Label.String = 'FRET Efficiency %';
    axis square;
    axis off;
    hold on;
    set(h, 'AlphaData', A);
    hold off

        saveas(fig5, [user.outputPath filesep user.saveName '_5.1 FRET Intensity Overlay' user.outputType])
    else
end


%% Figure 3. FRET Images of FA, Cytosol and Background (fixed min/max)
if flags.flagFRETimage || flags.flagOneFRETimage
    fig6 = figure("name", "FAs");
    imagesc(FRET_image_FA);
    colormap(BlackJet);
    title(Title_FA);
    c = colorbar;
    caxis([0 user.FRET_ColorbarLim]);
    c.Label.String = 'FRET Efficiency %';
    axis square;
    axis off;
    hold off;

    saveas(fig6, [user.outputPath filesep user.saveName '_5.2 FRET FA Image Fixed Range' user.outputType])
else
end

if flags.CytosolBackground == 1

    fig_cyto = figure("name", "Cytosol");
    imagesc(FRET_image_cyto);
    colormap(BlackJet);
    title(Title_cyto);
    c = colorbar;
    caxis([0 user.FRET_ColorbarLim]);
    c.Label.String = 'FRET Efficiency %';
    axis square;
    axis off;
    hold off;

    fig_bg = figure;
    imagesc(FRET_image_bg);
    colormap(BlackJet);
    title(Title_bg);
    c = colorbar;
    caxis([0 user.FRET_ColorbarLim]);
    c.Label.String = 'FRET Efficiency %';
    axis square;
    axis off;
    hold off;

else
end

if flags.CytosolBackground == 1
    saveas(fig_cyto, [user.outputPath filesep user.saveName '_2b. FRET Cytosol Image Fixed Range' user.outputType])
    saveas(fig_bg, [user.outputPath filesep user.saveName '_2c. FRET Background Image Fixed Range' user.outputType])
else
end


%% 11. Focal Adhesion Analysis NEW
if flags.FAMoreAnalysis

    % calculates the size of the pixel in micrometers
    fretPixelSizeMicrometer = user.Image_Size/user.Image_Resolution; 
    % get pixel area in micrometers²
    pixelAreaUm2 = fretPixelSizeMicrometer^2; %(should be in µm if defined as such (see inputs above)
    
    
    % switch for donor-acceptor as data source
    if data.flags.AcceptorData.switch
        % use acceptor channel for all intensity related metrics
        intensityImage = acceptor_channel;
    else
        % use donor (<Int>) channel for all intensity related metrics
        intensityImage = Int;
    end
    
    
    %%% Cell related metrics
    
    
    % compute cell area (in pixels and µm²)
    CellPixArea = sum(CellMask(:) == 1);
    cellAreaUm2 = pixelAreaUm2 * CellPixArea;
    % compute cell intensity
    totalCellIntensity = sum(intensityImage(CellMask == 1));
    % compute cytosol intensity (before removing low intensity background)
    cytosolIntensity = sum(intensityImage(CytosolMaskbefore ==1));
    
    
    %%% FA: Area
    
    
    % get area in pixels for each FA
    tmp = regionprops(faList, 'Area');
    data.fretAnalysis.AreaInPixels = [tmp.Area]';
    % calculates areas in micrometers for each FA (important for comparing different cells)
    data.fretAnalysis.AreaInMicroMeterSq = data.fretAnalysis.AreaInPixels * pixelAreaUm2 ;
    
    
    %%% FA: Intensity
    
    % compute overall FA intensity (integration over all FAs)
    data.fretAnalysis.totalFAintensity = sum(intensityImage(FAmask == 1));
    % calculates total intensity for each FA
    [data.fretAnalysis.TotalIntensity] = TotalValuePerArea(intensityImage, faList);
    data.fretAnalysis.TotalIntensity = data.fretAnalysis.TotalIntensity';
    % Calculates average intensity per FA & normalize
    data.fretAnalysis.AverageIntensity = AveragePerRegion(intensityImage, faList);
    data.fretAnalysis.AverageIntensityNorm = data.fretAnalysis.AverageIntensity / max(data.fretAnalysis.AverageIntensity, [], "all" );
    % calculates intensity per unit area (important when comparing different cells)
    data.fretAnalysis.IntensityPerMicroMeterSq = data.fretAnalysis.TotalIntensity ./ data.fretAnalysis.AreaInMicroMeterSq;
    % average FRET of all FAs combines
    % compute overall FA intensity (integration over all FAs)
    data.meanOfTotalFAintensity = mean(intensityImage(FAmask == 1), 'all');
    
    %%% FA: FRET
    
    numFAs = numel(faList.PixelIdxList);
    meanFRET = zeros(numFAs, 1);
    
    for i = 1:numFAs
        pixels = faList.PixelIdxList{i};
        values = FRET_image_FA(pixels);
        % Exclude NaN values
        values = values(~isnan(values));
        if ~isempty(values)
            meanFRET(i) = mean(values);
        else
            meanFRET(i) = NaN; % If all values are NaN, set mean to NaN
        end
    end
    
    statsMeanFRET = struct('MeanFRET', num2cell(meanFRET));
    meanFRETArray = [statsMeanFRET.MeanFRET];
    meanFRETList = meanFRETArray(:);
    data.fretAnalysis.AverageFRETperFA = meanFRET;
    
    %IF YOU WANT THE PREVIOUS VERSION, ACTIVATE THE FOLLOWING LINE AND
    %DEACTIVATE THE PREVIOUS LINE
    
    % assign Average FRET Index to each FAs.
    %data.fretAnalysis.AverageFRETperFA = AveragePerRegion(FRET_image_FA, faList);
    % normalize FRET data
    data.fretAnalysis.AverageFRETperFANorm = data.fretAnalysis.AverageFRETperFA / max(data.fretAnalysis.AverageFRETperFA, [], "all" );
else
end

%% COLOR AND INDEX OF FOCAL ADHESIONS

if flags.flagShowFaList
    figure("Name", 'Indexed FAs')
    PlotRegion(faList, faList.NumObjects, faList.PixelIdxList, faList.ImageSize, 'FA List (#', flags.flagTesting, 1, 1);
else
end

%Extra functions to find the index of FA by bringing mouse cursor on it.
%First select the Data Tips on top-right of the figure, then go to the FA
%and click on it!

if flags.flagFAindex == 1
    ShowFaIndex(faList)
else
end


%% FA Area vs Average and Total Intensity

if flags.flagShowFAAreaInt
    fig80 = figure("Name", 'Area vs Intensity');
    t = tiledlayout(1,2);

    % PANEL 1
    ax1 = nexttile;
    x = data.fretAnalysis.AreaInMicroMeterSq;
    y = data.fretAnalysis.TotalIntensity;

    scatter(x, y)
    scatter(x, y, 50, 'filled');
    hold on
    ScatterColors(x, y)
    hold on

    PlotLinReg(x, y, ':', 2, color.orangeTransp)
    corrA = corr(x, y);

    title(ax1, ['PCC = ' num2str(corrA)])
    ylabel(ax1, 'Total Intensity per FA');
    xlabel(ax1,'FA Area (µm²)');
    hold off

    % PANEL 2
    ax2 = nexttile;
    y = data.fretAnalysis.IntensityPerMicroMeterSq;

    scatter(x, y)
    scatter(x, y, 50, 'filled');
    hold on
    ScatterColors(x, y)
    hold on

    PlotLinReg(x, y, ':', 2, color.orangeTransp)
    corrB = corr(x, y);

    title(ax2, ['PCC = ' num2str(corrB)])
    ylabel(ax2, 'Average Intensity per FA');
    xlabel(ax2,'FA area (µm²)');
    hold off

        saveas(fig80, [user.outputPath filesep user.saveName '_4. Area vs Intensity' user.outputType])
    else
end

%% FRET vs Average and Total Intensity

if flags.flagShowFAFretInt
    fig90 = figure("Name", 'Intensity Plots');
    t = tiledlayout(1,3);

    % PANEL 1 Average Intensity vs FRET

    ax1 = nexttile;

    x = data.fretAnalysis.IntensityPerMicroMeterSq;
    y = data.fretAnalysis.AverageFRETperFA;
    scatter(x, y)
    scatter(x, y, 50, 'filled');
    hold on
    ScatterColors(x, y)
    hold on
    PlotLinReg(x, y, ':', 2, color.orangeTransp)
    corrC = corr(x, y);
    title(ax1, ['PCC = ' num2str(corrC)])
    ylabel(ax1, 'FRET Index');
    xlabel(ax1,'Average FA Intensity');
    hold off

    % PANEL 1 Total Intensity vs FRET

    ax2 = nexttile;

    x = data.fretAnalysis.TotalIntensity;
    y = data.fretAnalysis.AverageFRETperFA;
    scatter(x, y)
    scatter(x, y, 50, 'filled');
    hold on
    ScatterColors(x, y)
    hold on
    PlotLinReg(x, y, ':', 2, color.orangeTransp)
    corrC = corr(x, y);
    title(ax2, ['PCC = ' num2str(corrC)])
    ylabel(ax2, 'FRET Index');
    xlabel(ax2,'Total FA Intensity');
    hold off

    % PANEL 2 Area vs FRET

    ax3 = nexttile;
    x = data.fretAnalysis.AreaInMicroMeterSq;
    y = data.fretAnalysis.AverageFRETperFA;
    scatter(x, y)
    scatter(x, y, 50, 'filled');
    hold on
    ScatterColors(x, y)
    hold on
    PlotLinReg(x, y, ':', 2, color.orangeTransp)
    corrD = corr(x, y);
    title(ax3, ['PCC = ' num2str(corrD)])
    ylabel(ax3, 'FRET Index');
    xlabel(ax3,'FA Size (µm²)');
    hold off

        saveas(fig90, [user.outputPath filesep user.saveName '_5.FRET vs Intensity' user.outputType])
else
end

%% Line Analysis
% Cell shape metrics - centroid | mask | area (in pixels)


if flags.LineAnalysis % (decide to run this part or not)

    % toggle for generating the figures
    flags.ShowLineAnalysis = 0; % 0: skip, 1: show figures illustrating line analysis and orientation. 

    % ---------------------------------------------------------------------
    % Testing/validation: only keep region selected FAs
    flags.line.test = 0; % flags.line.test = 1;
    faLineList = faList;

    % ------------- Get shape metrics of FAs
    FA.shape.Circularity = regionprops(faLineList, "Circularity");
    FA.shape.Eccentricity = regionprops(faLineList, "Eccentricity");
    FA.shape.MajorAxisLengthPix = regionprops(faLineList, "MajorAxisLength");
    FA.shape.MinorAxisLengthPix = regionprops(faLineList, "MinorAxisLength");
    FA.shape.MajorAxisLengthUm = [FA.shape.MajorAxisLengthPix.MajorAxisLength]' .* fretPixelSizeMicrometer; % convert to actual length
    FA.shape.MinorAxisLengthUm = [FA.shape.MinorAxisLengthPix.MinorAxisLength]' .* fretPixelSizeMicrometer; % convert to actual length

    addpath('functions\LineAnalysis');
    IntMasked = double(Int);
    IntMasked(mask == 0) = NaN;
    IntMaskedResized = ResizeThis(IntMasked, data.logicalSize); 
    IntResized = ResizeThis(Int, data.logicalSize);
    MaskResized = ResizeThis(mask, data.logicalSize);

    FRETResized = ResizeThis(FRET_image_FA, data.logicalSize);
    falist_resized = bwconncomp(MaskResized);
    faLineList = falist_resized; 

    [CELL.centroid, CELL.mask, CELL.maskAreaPixels] = CentroidAreaOfCellFinder(IntResized);
    [FA.Perimeter] = PerimeterFA(CELL, faLineList, user.outputPath);

    % In this step, we do not define the orientation yet, we make a list of pixels of all perpendicular lines, but not from proximal to distal direction. 
    user.useCellCentroid = 1;
    [FA.intra.lineFaIdxList, FA.shape.centroidXY] = intraFAanalysisMain(faLineList, user, CELL);
    % VALIDATION if the orientation of FAs is correct, purely for visualization
    intraFAanalysisPlot01(flags.ShowLineAnalysis, faLineList, FA.intra.lineFaIdxList, FA.shape.centroidXY)
    % In this step, we finally include the orientation and make a list of pixels (proximal to distal) on each perpendicular line to the major axis. 
    [FA.intra.lineFaIdxListALL, FA.shape.centroidXYALL] = intraFAanalysisMain(faLineList, user, CELL);

    % ------------- Orientation of FAs
    [FA.orientation.vectorFromCenter, FA.orientation.polarHist.nBins, ...
        FA.orientation.polarHist.isBinNumber, FA.orientation.polarHist.binEdges, ...
        FA.Coordinates, FA.orientation.angleDegFA2Centroid] ...
        = OrientationFAtoCentroid(faList, FA.intra.lineFaIdxListALL, FA.shape.centroidXYALL);
    % VALIDATION 1: verify the orientation of FAs + centroids
    OrientationFAtoCentroid_figure1(CELL, flags.ShowLineAnalysis, faLineList, FA.shape.centroidXY, FA.Coordinates.FaCentroids, FA.Coordinates.ProximalVec, FA.Coordinates.DistalVec)
    % VALIDATION 2: verify the angle of FAs wrt centroid
    OrientationFAtoCentroid_figure2(flags.ShowLineAnalysis, FA.orientation.polarHist.binEdges, faLineList, FA.orientation.polarHist.isBinNumber); % check FA orientation bins
    
    % Orientation vs Distance
    FA.Coordinates.vectorFromCenterNormPix = OrientationDistance(flags.ShowLineAnalysis, ...
        FA.orientation.vectorFromCenter, FA.orientation.polarHist.nBins, ...
        FA.orientation.polarHist.isBinNumber, FA.orientation.polarHist.binEdges);

    user.line.interpolationSampleSize = 50;
    % ------------- LINE ANALYSIS : FRET --------
    [FA.intra.FRET, FA.intra.FRETaverage] = intraFAanalysisExtract(faLineList, FA.intra.lineFaIdxList, FRETResized);
    [FA.intra.interp.FRET, FA.intra.interp.FRETaverage, FA.intra.interp.normFRET, FA.intra.interp.normFRETaverage] = ...
        interpolateLineAnalysis('FRET', FA.intra.FRETaverage, user.line.interpolationSampleSize, user);
    [FA.intra.FRETpadded] = LineAnalysisPadding (FA.intra.FRETaverage, FA.intra.interp.FRET, user.line.interpolationSampleSize);
    % ------------- LINE ANALYSIS : INTENSITY ------
    [FA.intra.intensity, FA.intra.intensityAveraged] = intraFAanalysisExtract(faLineList, FA.intra.lineFaIdxList, IntMaskedResized);
    [FA.intra.interp.intensity, FA.intra.interp.intensityaverage, FA.intra.interp.normIntensity, FA.intra.interp.normIntensityaverage] = ...
        interpolateLineAnalysis('Intensity',FA.intra.intensityAveraged, user.line.interpolationSampleSize, user);
    [FA.intra.Intensitypadded] = LineAnalysisPadding (FA.intra.intensityAveraged, FA.intra.interp.intensity, user.line.interpolationSampleSize);
    
end

%% Extra metrics

% need to position this somewhere else
MaskFromFRET = ~isnan(FRET_image_FA); % why again? because additonal regions got removed (e.g. value too low) -> this is to make it match again
avFRET = str2double(FRETav);

% vars to keep
%outputMetrics = data.fretAnalysis;
outputMetrics.Info = user; % what parameters were used

outputMetrics.FRETimageFA = FRET_image_FA;
outputMetrics.Mask = MaskFromFRET;
outputMetrics.FRETav = avFRET;

outputMetrics.IntD = Int;
outputMetrics.IntA = acceptor_channel;

outputMetrics.dataPath = data.dataPath;

%outputMetrics.meanOfTotalFAintensity = data.meanOfTotalFAintensity;
%outputMetrics.FRETav_bg = FRETav_bg;
%outputMetrics.FRETav_cyto = FRETav_cyto;
%outputMetrics.FRETindex = FRET_Index_Data;
%outputMetrics.FRETlist = FRET_list;
%outputMetrics.Mask = FAmask;

%outputMetrics.cellAreaUm2 = cellAreaUm2;
%outputMetrics.CellPixArea = CellPixArea;
%outputMetrics.pixelAreaUm2 = pixelAreaUm2;
%outputMetrics.RatioFACellArea = sum(outputMetrics.AreaInMicroMeterSq)/outputMetrics.cellAreaUm2*100;
%outputMetrics.cellIntensity = totalCellIntensity;
%outputMetrics.cytosolIntensity = cytosolIntensity;
%outputMetrics.AreaInPixelsFA = data.fretAnalysis.AreaInPixels;

if flags.flagPhasorHistogram
    outputMetrics.Histogram_BinCenters = binCenters;
    outputMetrics.Histogram_Counts = FRET_counts;
    outputMetrics.Histogram_NormalizedCounts = FRET_counts_norm;
end

end
