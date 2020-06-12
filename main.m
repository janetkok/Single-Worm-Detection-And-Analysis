% Kok Yong En
% hfyyk1
% 10346769
% Version 2
% 30/5/2020

close all;
clear;
clc;

choice = menu('Load from','specific images','img numbered 1:7');
switch choice
    case 1
        %% load from 1 or more specific images
        try_again = true;
        while try_again %user must select at least 1 image
            [FileNames,PathName] = uigetfile('*.png', 'Chose images to load:','MultiSelect','on');
            nfiles = length(FileNames);
            if nfiles == 1%if there is no image selected bcuz if 1 image is selected then it will be char instead
                nfiles = 0;
            end
            if nfiles ~= 0; try_again = false; end
            if try_again
                warning('Please select at least 1 image')
            end
        end
        
        if isa(FileNames,'char')%if 1 image is selected
            images = strings([1]);
            images(1) = fullfile(PathName,FileNames);
        else %if multiple images are selected
            images = strings([1,nfiles]);
            for i = 1:nfiles
                images(i) = fullfile(PathName,FileNames{i});
            end
        end
        
    case 2
        %% Array of 1:7 images from img folder
        num = string([1:7]);
        images = strcat('img/',num,'.png');
    otherwise
        disp('closing program');
        return
end

%% Loop through all images
for image_url = images
    % Read image
    im = imread(image_url);
    
    %% Contrast stretching,Sharpening,Grayscale conversion
    st = imadjust(im,stretchlim(im));
    sharp = imsharpen(st,'Radius',1.5,'Amount',2);
    im_gray = rgb2gray(sharp);
    
    %% Edge detection and subtract darker background
    % detect worm edge using log operator
    [~,threshold] = edge(im_gray,'log');
    fudgeFactor = 2.2;
    log = edge(im_gray,'log',threshold * fudgeFactor);
    
    % detect worm edge using Prewitt operator
    [~,threshold] = edge(im_gray,'Prewitt');
    fudgeFactor = 0.9;
    prew = edge(im_gray,'Prewitt',threshold * fudgeFactor);
    prew = bwareaopen(prew,3);
    
    %detect unnecessary lines of dark background and petri dish
    background = edge(im_gray,'Canny',0.55);
    
    %subtract darker background for log edge detection
    se1=strel('disk',1);
    backgroundD1 = imdilate(background,se1);
    logClean = log&(~backgroundD1);
    
    %subtract darker background for prewitt edge detection
    se4=strel('disk',4);
    backgroundD4 = imdilate(background,se4);
    prewClean = prew&(~backgroundD4);
    
    %background darkerness of image
    green_chan = im(:,:,2);
    blue_chan = im(:,:,3);
    gb = (green_chan/2 + blue_chan/2);
    % Otsu thresholding on gb
    t = graythresh(gb);
    middleground = imbinarize(im_gray, t);
    middlegroundClean = bwareaopen(~middleground,10);
    se4 = strel('disk',4);
    middlegroundCClosed = imclose(middlegroundClean,se4);
    middlegroundCCC = bwareaopen(middlegroundCClosed,130);
    %subtract middleground
    prewClean2 = prewClean&(~middlegroundCCC);
    
    %combine both cleaned binary gradient masks
    BW = logClean|prewClean2;
    
    %% Filling holes of worm(BW) with cleaned hole fillers of dilated worm(BWdil)
    % get clean dilated holes fillers by filling holes of dilated worm
    % BWdil mask. then use clean dilated hole fillers to fill holes of BW
    % mask(not dilated)
    
    %perpendicular linear structuring element
    se90 = strel('line',2,90);
    se0 = strel('line',2,0);
    %first trial of dilate worm mask to fill gaps
    BWdil = imdilate(BW,[se90 se0]);
    %first trial of filling holes of dilated worm mask
    BWdfill = imfill(BWdil,'holes');
    %get hole fillers
    holes = BWdfill - BWdil;
    %clean hole fillers
    holesClean = bwareaopen(holes,2);
    %dilate cleaned hole fillers
    se1 = strel('disk',1);
    holesCDil = imdilate(holesClean,se1);
    
    %filling BW mask(not dilated) with clean dilated hole fillers
    BWfill2 = holesCDil|BW;
    
    
    %% Noise removal and fill back holes loss during cleaning of hole fillers
    % Retain only 1 largest blob to remove noise
    BWclean = bwareafilt(BWfill2, 1,'largest');
    
    %dilate hole fillers
    holesDil = imdilate(holes,se1);
    %combine dilated hole fillers with segmented blob
    BWcleanFill = holesDil|BWclean;
    
    %Retain 1 largest blob
    BWC = bwareafilt(BWcleanFill, 1,'largest');
    
    %% Fill gaps and holes, smoothing
    %dilate worm lines to fill gaps
    BWCdil = imdilate(BWC,[se90 se0]);
    se6=strel('disk',6);
    %closing and opening to fill holes
    BWCclosed = imclose(BWCdil,se6);
    se7=strel('disk',7);
    BWCopened = imopen(BWCclosed,se7);
    %smoothing the worm
    se2 = strel('disk',2);
    worm = imerode(BWCopened,se2);
    
    %% Area of worm
    area = bwarea(worm);
    
    %% Boundary
    % Trace boundary around worm
    edges = edge(worm ,'canny');
    se = strel('square', 2);
    highlights = imdilate(edges, se);
    
    %% Skeleton, length, tortuosity
    % skeletonize
    skelImageF = bwskel(worm); % First guess with short spurs.
    %take 1/3 of first guess spine length as min spine length
    minBranchLength = round(sum(skelImageF(:))/3);
    % prune all spurs shorter than minBranchLength
    skelImage = bwskel(worm, 'MinBranchLength', minBranchLength);
    
    % spine length.
    spineLength = sum(skelImage(:));
    
    % endpoints of spine
    endpointImage = bwmorph(skelImage, 'endpoints');
    [endR, endC] = find(endpointImage);
    % straight line distance of endpoints
    straightLineDistance = sqrt((endC(2) - endC(1))^2 + (endR(2) - endR(1))^2);
    % tortuosity of worm
    tortuosity = spineLength / straightLineDistance;
    
    %% mean and middle radius, mean diameter,volume
    % Get the Euclidean Distance Transform.
    edtImage = bwdist(~worm);
    % Measure mean radius by looking along the skeleton of the distance transform.
    meanRadius = mean(edtImage(skelImage));
    % mean diameter
    meanDiameter = 2 * meanRadius;
    
    %coordinates of skeleton
    skelCo = bwtraceboundary(skelImage,[endR(1) endC(1)],'S');
    skelR = skelCo(:,1); %skeleton rows
    skelC = skelCo(:,2);%skeleton columns
    %coordinates of boundaries
    boundariesCo = bwboundaries(worm);
    boundariesCo = boundariesCo{1};
    boundR = boundariesCo(:, 1); % boundary rows
    boundC = boundariesCo(:, 2);% boundary columns
    
    % middle index of skeleton
    midIndex = ceil(spineLength/2);
    % Find distances from mid of skeleton coord(x, y) to all other boundary coordinates.
    distances = sqrt((boundC - skelC(midIndex)).^2 + (boundR - skelR(midIndex)).^2);
    % middle radius(shortest distance from boundary coord) and index of closest boundary coord
    [midRadius, indexOfMin] = min(distances);
    % store the closest boundary coordinates
    midC = boundC(indexOfMin);
    midR = boundR(indexOfMin);
    
    %radii coordinates matrix initialization
    radBoundC = zeros(spineLength,1);
    radBoundR = zeros(spineLength,1);
    %length of radii matrix
    lenRadius  = ceil(spineLength/10);
    %radii matrix initialization
    radii = zeros(lenRadius,1);
    
    %radii at interval of 10 pixels along skeleton
    for j = 1:lenRadius
        n = 10*j - 9; %skeleton index at intervals of 10pixels
        % Find distances from skeleton coordinate(x, y) to all other boundary coordinates.
        distances = sqrt((boundC - skelC(n)).^2 + (boundR - skelR(n)).^2);
        % Find the min radius and index of closest boundary coordinate
        [minRadius, indexOfMin] = min(distances);
        % store radius of closest distance from boundary coordinate
        radii(j) = minRadius;
        % store the closest boundary coordinates
        radBoundC(j) = boundC(indexOfMin);
        radBoundR(j) = boundR(indexOfMin);
    end
    
    %volume
    volume = 0.0;
    for i=1:(length(radii)-1)
        r1 = radii(i);
        r2 = radii(i+1);
        volume = volume + ((pi*10/3)*(r1^2 + r1*r2 + r2^2));
    end
    
    %% Display results
    figure('name',"Output for "+ image_url);
    
    % contrast stretching, shrapening, grayscale conversion
    subplot(3,3,1);
    imshow(im_gray)
    title('Contrast stretching, Shrapening, Grayscale conversion','FontSize', 9,'Interpreter', 'None');
    
    % edge detection and subtract darker background
    subplot(3,3,2);
    imshow(BW)
    title('Edge detection and subtract darker background','FontSize', 9,'Interpreter', 'None');
    
    % Filling holes of worm(BW) with cleaned hole fillers of dilated worm(BWdil)
    subplot(3,3,3);
    imshow(BWfill2)
    title('Filling holes of worm with cleaned hole fillers of dilated worm','FontSize', 9,'Interpreter', 'None');
    
    % Noise removal and fill part loss during cleaning of hole fillers
    subplot(3,3,4);
    imshow(BWC)
    title('Noise removal and fill part loss during cleaning of hole fillers','FontSize', 9,'Interpreter', 'None');
    
    % Fill gaps and holes, smoothing
    subplot(3,3,5);
    imshow(worm);
    title('Fill gaps and holes, smoothing','FontSize', 9,'Interpreter', 'None');
    
    % overlay segmented worm mask, calculate area
    subplot(3,3,6);
    wormMask = labeloverlay(im,worm,'Colormap','cool','Transparency',0.7);
    imshow(wormMask)
    caption = sprintf('Area: %.2f',area);
    title(caption,'FontSize', 10,'Interpreter', 'None');
    
    % trace boundary, compute skeleton length, tortuosity
    subplot(3,3,7);
    highlightedWorm = imoverlay(im,highlights,'cyan');% boundary
    boundarySkel = imoverlay(highlightedWorm,skelImage,'green');%skeleton
    imshow(boundarySkel);
    axis('on', 'image');
    hold on;
    dummyCyan = plot(boundC(1),boundR(1), 'Color', 'c'); %dummy plot for cyan line
    dummyGreen = plot(endC(1),endR(1), 'Color', 'g'); %dummy plot for green line
    sld = plot(endC, endR, 'r-', 'LineWidth', 2); %plot straightLineDistance of endpoints
    title('trace boundary, compute skeleton length, tortuosity', 'FontSize', 9, 'Interpreter', 'None');
    spineCap = sprintf('Spine Length:%d', spineLength);
    tortuosityCap = sprintf('Tortuosity:%.2f',tortuosity);
    legend([dummyCyan,dummyGreen,sld],{'boundary',spineCap,tortuosityCap},'Location','northwest','FontSize',8)
    
    % euclidean distance transform image,mean radius, mean diameter
    subplot(3,3,8);
    imshow(edtImage,[]); %euclidean distance transform image
    caption = sprintf('Euclidean Distance Transform. Mean radius: %.2f Mean diameter: %.2f',meanRadius,meanDiameter);
    title(caption, 'FontSize', 8, 'Interpreter', 'None');
    
    % middle radius,volume
    subplot(3,3,9);
    highlightedWorm = imoverlay(im,highlights,'cyan');% boundary
    boundarySkel = imoverlay(highlightedWorm,skelImage,'yellow');%skeleton
    imshow(boundarySkel);
    axis('on', 'image');
    hold on;
    for j = 1:lenRadius %plot radius at interval of 10 pixels along skeleton
        n = 10*j - 9;
        rad = plot([radBoundC(j),skelC(n)], [radBoundR(j),skelR(n)], 'g-', 'LineWidth', 1);
    end
    mid = plot([midC,skelC(midIndex)], [midR,skelR(midIndex)], 'b-', 'LineWidth', 2); %plot middle radius
    caption = sprintf('Middle radius:%.2f Volume:%.2f', midRadius,volume);
    title(caption, 'FontSize', 9, 'Interpreter', 'None');
    legend([rad mid],{'Radius at interval of 10','Middle radius'},'Location','southwest','FontSize',7)
    
    fprintf('processed %s\n',image_url);
end
