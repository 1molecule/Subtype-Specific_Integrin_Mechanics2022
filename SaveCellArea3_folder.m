% Detect and save the cell area (binary) using Image Segmentation
% optimized for RICM images
% MJ 2019-2021
% It analyzse all images in the assigned folder or all subfolders.
% 2022. 4: global thresholding result added. zero image added. 
% 2022.5.1.: adaptive thresholding adding...

clear all
close all

%% Settings
fAllsubfolders = true;      % "false" for the selected sigle folder. 1/1.04
chan = 2;                   %my RICM is in the 2nd channel.
fEdgeBased = true;          %area dection method based on edge dection
fGlobalThBased = true;      %area dection method based on global thres.
fAdaptiveThBased = 0;       %area dection method based on adaptive thres.
fSaveBlank = 0;             %Save the all zero image
fSaveTheChan = true;        %Save the auto-contrasted one-chan image
extension = 'tif';
savePathTag = '_mask';      %make subfolder to save the result or []
fAutoContrast = false;      %Use "imadjust" before cell region detection
MedianFilterMat = [0 0];    %0 not to use; Can I remove the shadow of spherial cell? No.

%Parameters for the segmentaion-based method.
fCheck_seg = false;
dialation = 4;              %Dilation factor, usually 5
fNoBorder = fAllsubfolders; %Remove regions touching the edge of image.
%fNoBorder = false;
fLargest = false;           % Select the largest region only
method_edge = 'Canny';      %'sobel' or 'Canny' etc.
%fudgeFactorVec = [0.5, 0.8, 1, 1.5, 2, 2.5, 3, 4, 5];%w/ Canny for RICM
fudgeFactorVec = [ 0.8, 1, 1.2, 1.5, 1.8, 2, 2.5, 3, 4];%w/ Canny for RICM 2
%fudgeFactorVec = [1.8, 2, 2.2, 2.5, 2.8, 3, 4, 5, 6];%w/ Canny for RICM 2
%fudgeFactorVec = [1, 1.2, 1.5, 1.8, 2, 2.2, 2.5, 2.8, 3, 4, 5]; % test
%fudgeFactorVec = [0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9]; %w/ sobel for RICM
%fudgeFactorVec = [0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7]; % DIC
%FudgeFactor = 4; %% DIC 0.2 / RICM 1 ??

%Parameters for the thresholding-based method.
fCheck_thr = false;
nLargestVec = [1, 2, 3];   %choose n largest region.


%% Find the path & subfolders (options: ending with '_each')
path0 = uigetdir('C:\Aworkspace');

if fAllsubfolders
    files    = dir(path0);
    names    = {files.name};
    dirFlags = [files.isdir] & ~strcmp(names, '.') & ~strcmp(names, '..');
    subfolders = names(dirFlags);
    numFolder = size(subfolders,2);
else
    numFolder = 1;
end

%% Anaylize the subfolders
for cnt_subfoler = 1:numFolder
    
    if fAllsubfolders
        parentFolder = path0;
        path2 = subfolders(cnt_subfoler);
        thisFolder = path2{1}
        path = [path0 filesep thisFolder];
    else
        path = path0;
        [parentFolder,thisFolder,~]=fileparts(path);
    end
    

    %% Find all image stacks in the folder
    path_save = [path0 filesep thisFolder savePathTag];
    status = mkdir(path_save);
    cd(path);
    %matfiles = dir(['**/*.' extension]);
    matfiles = dir(['*.' extension]);
    %matfiles = dir(fullfile(thisPath, '*.tif'));
    nfiles = length(matfiles);
    
    for i = 1:nfiles
        
        filename = matfiles(i).name
        
        if contains(filename,'region')
            continue
        end
        
        if extension == 'nd2'
            data1 = bfopen(filename);
            [size1, size2] = size(data1{1, 1}{1, 1});
            img = data1{1, 1}{chan, 1};
        else
            img = imread(filename, chan);    % for TIFF
        end
        
        %filenamehead1 = filename(1:end-4);
        [dummy,filenamehead1,ext] = fileparts(filename);
        
       %% image analysis       
        if fAutoContrast
            img = imadjust(img);
        end
        
        if MedianFilterMat(1)
           img =  medfilt2(img,MedianFilterMat);
        end
        
        if fEdgeBased
            % 1. Segmetaion-based, for multiple fudge factors
            for j=1:size(fudgeFactorVec,2)
                % analysis & save
                fudgeFactor = fudgeFactorVec(j);
                
                regions = mjf_cell_regions(img, method_edge, fudgeFactor, dialation, fNoBorder, fCheck_seg);
                filename4save = [path_save filesep filenamehead1 '_region_ff' ...
                    num2str(fudgeFactor) '.tif'];
                if fLargest
                    region1 = bwareafilt(regions,1);
                    imwrite(region1, filename4save);
                else
                    imwrite(regions, filename4save);
                end
            end
        end
        
        if fGlobalThBased
            % 2. theshold-based *n largest only...
            for j=1:size(nLargestVec,2)
                nLargest = nLargestVec(j);
                regions = mjf_cell_regions_theshold(img, nLargest, fCheck_thr);
                filename4save = [path_save filesep filenamehead1 '_region_thes' ...
                    num2str(nLargest) '.tif'];
                imwrite(regions, filename4save);
            end
            
        end
        
        if fAdaptiveThBased
            % 3. adaptive thresholding based method
            regions = mjf_cell_regionsr_aTh(img);
            filename4save = [path_save filesep filenamehead1 '_region_aTh.tif'];
            imwrite(regions, filename4save);
        end
        
        if fSaveBlank
            % save the empty image.
            non = zeros(size(img), 'logical');
            filename4save = [path_save filesep filenamehead1 '_region_0.tif'];
            imwrite(non, filename4save);
        end
        
        % save the original image (adjusted) for comparision
        if fSaveTheChan
            if ~fAutoContrast
                Ia = imadjust(img);
            else
                Ia = img;
            end
            imwrite(Ia, [path_save filesep filenamehead1 '_AC_chan' num2str(chan) '.tif']);
        end
        %fclose(data1);
    end
    
end

cd(parentFolder);
% clear all
% fclose all


%%     

function BW = mjf_cell_regionsr_aTh(img)
sensitivity = 0.56;         % use 0.55-"0.56" for RICM; 0.5 by MATLAB default.
nbhd1 = [7 199];            % [3 199] 2-element vector of positive odd integers
ThSmooth = 32;              % 64 smooth the threshold matrix to minimize the effect of "the regions"
sizeFilter = [9 inf];      % pixel; size filter for binary image
MedianFilterMat = [3 3];    % [5 5]
fFilling = false;           % It is more consitent with TH method not to fill.
nCheck = 0; check=1;        % 7; 0 not to use, show the result of each step
fNoBorder = true;

if nCheck; subplot(1, nCheck, check); imshow (img, []); title('img'); colormap jet; check= check+1; end

img =  medfilt2(img,MedianFilterMat);
if nCheck; subplot(1, nCheck, check); imshow (img, []); title('Median'); check= check+1; end

T1 = adaptthresh(img, sensitivity,'ForegroundPolarity','dark', 'NeighborhoodSize', nbhd1);
if ThSmooth
    T1 = imgaussfilt(T1, ThSmooth);
end

BW = imbinarize(img,T1);
BW = imcomplement(BW);  % invert the image, we use dark region of RICM
if nCheck; subplot(1, nCheck, check); imshow (BW, []); title('TH'); check= check+1; end

% BW =  medfilt2(BW,MedianFilterMat);
% if nCheck; subplot(1, nCheck, check); imshow (BW, []); title('Median'); check= check+1; end

dialation = 3;
se90 = strel('line', dialation, 90);
se0 = strel('line', dialation, 0);
BW = imdilate(BW, se90);
BW = imdilate(BW, se0);

seD = strel('diamond',1);
BW = imerode(BW,seD);
%BW = imerode(BW,seD);
%imerode is not good for transient adhesion with many small compenents.

BW = bwareafilt(BW, sizeFilter); % size filter
if nCheck; subplot(1, nCheck, check); imshow (BW, []); title('sizeFilter'); check= check+1; end

if fNoBorder
    BW = imclearborder(BW, 4);
end

if fFilling
    BW = imfill(BW,'holes');
    if nCheck; subplot(1, nCheck, check); imshow (BW, []); title('imfill'); check= check+1; end
end

%BWoutline = bwperim(BWfinal);

end

function BW = mjf_cell_regions_theshold(I, nLargest, fCheck)
% [filename,path] = uigetfile('*.tif');
% cd(path);
% chan=2;
% I = imread(filename, chan);

if fCheck
figure, tiledlayout('flow'), nexttile;
imshow(I, []), title('original image');
end 

BW = imcomplement(imbinarize(I));
if fCheck; nexttile; imshow(BW); title('global'); end

%BWnobord = imclearborder(BW_erode1, 4);
BWnobord = imclearborder(BW, 4);
if fCheck; nexttile; imshow(BWnobord); title('BWnobord'); end

seD = strel('diamond',1);
BW_erode1 = imerode(BWnobord,seD);
if fCheck; nexttile; imshow(BW_erode1); title('erode1'); end

BW_erode2 = imerode(BW_erode1,seD);
if fCheck; nexttile; imshow(BW_erode2); title('erode2'); end

region1 = bwareafilt(BW_erode2,nLargest);
if fCheck; nexttile; imshow(region1); title('largest'); end

dialation = 3;
se90 = strel('line', dialation, 90);
se0 = strel('line', dialation, 0);
region1_pre = imdilate(region1, se90);
region1_dia = imdilate(region1_pre, se0);
if fCheck; nexttile; imshow(region1_dia); title('dialate'); end

filled = imfill(region1_dia, 'holes');
if fCheck; nexttile; imshow(filled); title('filled'); end

BW = filled;
%BWoutline = bwperim(BWfinal);
end

function BW = mjf_cell_regions(I, method_edge, fudgeFactor, dialation, fNoBorder, fCheck)
% MJ 171126 - check cell_boundary_seg to test the code
% Based on "Detecting a Cell Using Image Segmentation"
% https://www.mathworks.com/help/images/examples/detecting-a-cell-using-image-segmentation.html

%method_edge = 'Canny';      %'sobel' 'Prewitt' 'Roberts' 'log' 'zerocross' 'Canny' 'approxcanny'
%len = 5; %Dilation factor

%[regions, outlines] = mjf_cell_regions(I, method_edge, fudgeFactor, dialation, fNoBorder);
%Step 1: Read Image
%I = stack1(corner(3,2):corner(4,2),corner(1,1):corner(2,1),1)

if fCheck
    figure, tiledlayout('flow'), nexttile;
    imshow(I, []), title('original image');
end

%% Additonal filtering (Step 1.5)
I = medfilt2(I);

%% Step 2: Detect Entire Cell
[~, threshold] = edge(I, method_edge);
BWs = edge(I, method_edge, threshold*fudgeFactor);
if fCheck
%     for z=1:10
%         fudgeFactor = 0.4*z;
%         BWs = edge(I, method_edge, threshold*fudgeFactor);
%         nexttile, imshow(BWs), title(num2str(fudgeFactor));
%     end
    BWs = edge(I, method_edge, threshold*fudgeFactor);
    nexttile, imshow(BWs), title(num2str(fudgeFactor));
end

%% Step 3: Dilate the Image
%len = 5;
se90 = strel('line', dialation, 90);
se0 = strel('line', dialation, 0);
BWsdil_pre = imdilate(BWs, se90);
BWsdil = imdilate(BWsdil_pre, se0);
if fCheck
    nexttile, imshow(BWsdil_pre), title('Dilate90');
    nexttile, imshow(BWsdil), title('Dilate90and0');
end

%% Step 4: Fill Interior Gaps
BWdfill = imfill(BWsdil, 'holes');
if fCheck
    nexttile, imshow(BWdfill), title('image filled');
end

%% Step 5: Remove Connected Objects on Border
if fNoBorder
    BWnobord = imclearborder(BWdfill, 4);
else
    BWnobord = BWdfill;
end
if fCheck
    nexttile, imshow(BWnobord), title('border option');
end

%% Step 6: Smoothen the Object
seD = strel('diamond',1);
BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWfinal,seD);

%BWoutline = bwperim(BWfinal);
%regions = BWfinal;
%outlines = BWoutline;

if fCheck
    outlines = bwperim(BWfinal);    
    nexttile, imshow(outlines), title('smoothened result');
    %temp = input('fCheck is on. Enter to proceed.');
end

BW = BWfinal;

end