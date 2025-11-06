% currentFolder = 'C:\LCSD\P63 dataset 111821\P63 dataset 111821\p63 CLIN_0009\D2 overlays'
% img = imread("0_D1_R1_01_Overlay.tif");
% cp = cellpose(Model="cyto2");
% imgB= img(:,:,3);
% 
% labels = segmentCells2D(cp,imgB,ImageCellDiameter=7);
% imgOvrlay = labeloverlay(imgB,labels);
% imshow(imgOvrlay)

%   =======================================================================================
%   Copyright (C) 2013  Erlend Hodneland
%   Email: erlend.hodneland@biomed.uib.no 
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   =======================================================================================

clear all
close all
clc
% problems: CLIN_0008,D1, CLIN_0012,D3 (image size is  4x larger), Entire
% CLIN_0015. CLIN_0006 D1,D2,D3 directories missing, CLIN_0005 having
% problems too, CLIN_0004 has larger images too but does have a downsample
% directory.
AnalysisList = readtable('C:\LCSD\P63 dataset 111821\AnalysisBookForSegmentationCompareTemp.xlsx','Format','auto');
DataIdx1 = 1;
for K = 1:size(AnalysisList.Folder,1)

currentFolder = AnalysisList.Folder{K}
%currentFolder = 'C:\LCSD\P63 dataset 111821\P63 dataset 111821\P63 CLIN_0009\D1 overlays'
CurrDir= currentFolder(end-20:end-12); % Check CLIN_0012, D3 later for large image size
statsDir = currentFolder(1:end-12);
statsSubDir = currentFolder(end-10:end-9);
WriteDir = string('C:\LCSD\Results\')+CurrDir;
CurrDir = insertBefore(CurrDir,"_","\"); % Insert \ before underscore to display correctly in subtitle
StoreImages = 1;
StoreImagesTif = 1;
StoreImageSidebySide = 0;
StoreImageSidebySide1 = 1;
StoreForTraining = 0;
StoreForTraining1 = 0;
DisplayImg = 1;

% Is image lower resolution/ (940x780)
SF = 1;

cd(currentFolder)
%addpath(currentFolder)
ImageList = natsortfiles(ls('*.tif'));

% % % load the data
% % load ../data/nucleistain_3D.mat
% % imsegmB = imsegmB(:,:,7);
%%%SS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for imInx = 1:size(ImageList,1)

% Is image lower resolution/ (940x780)
SF = 1;   
cp = cellpose(Model="cyto");
%cp = cellpose(Model="cpsam");
%cp = cellpose(Model="cyto",UseEnsemble=true);

    toggle = 1;
%currentFolder = 'C:\Users\junsu\Dropbox\LCSD\P63 dataset 111821\P63 dataset 111821\p63 CLIN_0009\D1 overlays'
% cd(currentFolder)
%addpath(currentFolder)

theimage = strtrim(ImageList(imInx,:))
%theimage = 'D3_R1_01_Overlay.tif'; % This is only when you want to force reading an image, also force current_folder in line 32
if theimage == "_"
  ImgTag=theimage(1:end-12);
else
   ImgTag=theimage(1:end-11);
end
Wfile = theimage(1:end-12);
ImgTag = insertBefore(ImgTag,"_","\"); % Insert \ before underscore to display correctly in subtitle

Stat = [];

I = imread(theimage);
I2 = im2double(I);
imsegmB = I2(:,:,3); 

if DisplayImg == 1
  cellsegm.show(imsegmB,1);
  title({CurrDir+": Raw Image";ImgTag});axis off;
  figure; figure(123); imshow(imsegmB)
end
if StoreImages
    saveas(gcf, WriteDir+string('\')+string('\Orig\')+string(Wfile), 'pdf');
end
if StoreImagesTif
    imwrite(mat2gray(imsegmB), WriteDir+string('\')+string('OrigTif\')+string(Wfile)+string('.tif'));
end
%(101);
%imshow(I);

%cellbwCellseg = segmentCells2D(cp,imsegmB,ImageCellDiameter=9); %Default
%se = strel('disk', 10);
%tophatFiltered = imtophat(imsegmB,se);
%tophatFiltered = imadjust(tophatFiltered);
cellbwCellseg = segmentCells2D(cp,imsegmB,ImageCellDiameter=10); % Default
%cellbwCellseg = segmentCells2D(cp,imsegmB,ImageCellDiameter=10, tile=true, ...
%    TileOverlap=0.1);
cellbw1 = cellbwCellseg;
cellbw1stat = bwconncomp(cellbw1);
CC=imoverlay(imsegmB,cellbw1);
%imshow(labeloverlay(imsegmB,cellbw1))
% 
imageNames = theimage;
cellCounts1 = cellbw1stat.NumObjects;
if DisplayImg == 1
  cellsegm.show(CC,2);
  title({CurrDir+": Cell segmentation, BLUE: "+ cellbw1stat.NumObjects + " Cells";ImgTag})
  axis off;
end
if StoreImages
%    saveas(gcf, WriteDir+string('\')+string('\Segmented\')+string(Wfile), 'tif');
end
%figure;
%imshow(label2rgb(cellbwCellseg, 'jet', 'k', 'shuffle'));
%title('Segmented Cells');
overlayImage = imoverlay(imsegmB, cellbw1, 'red');
% Get the boundaries of each cell region
props = regionprops(cellbw1, 'PixelIdxList', 'BoundingBox', 'Perimeter');

% Get the perimeter of each cell region in cellbw1
boundaryPixels = bwperim(cellbw1); % Extract only the boundaries

% Thicken the boundary by a small amount (1-pixel thickening)
thickBoundaryPixels = imdilate(boundaryPixels, strel('disk', 1)); % Adjust the size of the disk if needed

% Define a bright red boundary color (higher intensity)
brightRed = [1, 0.2, 0.2];

% Overlay the thickened boundary on the original RGB image I2
boundaryOverlay = imoverlay(I2, thickBoundaryPixels, brightRed);

% Display the boundary overlay image
if DisplayImg == 1
    imshow(boundaryOverlay);
    title({CurrDir+", "+ ImgTag+": Cell Boundaries in Bright Red: "+ cellbw1stat.NumObjects + " Cells"});
end
% Save the displayed overlaid image
if StoreImages
    % Define the absolute path for the output folder
    outputFolder = fullfile('C:\LCSD\Results', regexprep(CurrDir, '[\\/:*?"<>|]', ''), 'Segmented'); 
    % Create the output folder if it does not exist
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    % Clean ImgTag by replacing special characters with underscores
    cleanImgTag = regexprep(Wfile, '[\\/:*?"<>|]', '');

    % Construct the fully qualified file path
    outputFileName = fullfile(outputFolder, strcat(cleanImgTag, '.tif'));

    % Display the output filename for verification
    disp(['Saving to: ' outputFileName]);

% Create a new figure to display the image and set the title
    figureHandle = figure;
    imshow(boundaryOverlay);
    title([string(CurrDir) + ", " + string(ImgTag) + ": Cell Boundaries in Bright Red: " + num2str(cellbw1stat.NumObjects) + " Cells"]);

    % Save the figure with the title as a .tif image
    exportgraphics(figureHandle, outputFileName, 'Resolution', 300);

    % Close the figure to avoid clutter
    close(figureHandle);
end

%Splitting cells using the watershed transform
%Follows the notes in https://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/
% Convert the labeled image to binary
binaryCellSeg = cellbwCellseg > 0;

D = -bwdist(~binaryCellSeg);
D(~binaryCellSeg) = Inf;
%imshow(D,[])
mask = imextendedmin(D,1);
%imshowpair(cellbwCellseg,mask,'blend')
D2 = imimposemin(D,mask);
%D3 = imhmin(D2,5); %10 is the height threshold for suppressing shallow minima
Ld2 = watershed(D2);

bw3 = cellbwCellseg;
bw3(Ld2 == 0) = 0;
cellbw3stat = bwconncomp(bw3);
%bw3_binaryMask = bw3>0;

if DisplayImg == 1
  cellsegm.show((bw3>0),3);
  title({CurrDir+": Segmentation (post cell split), BLUE: "+ cellbw3stat.NumObjects + " Cells";ImgTag})
  axis off;
end

bw3 = checknghbr(bw3);
cellbw3stat = bwconncomp(bw3);
bw3_binaryMask = bw3>0;
cellCounts2 = cellbw3stat.NumObjects;
if DisplayImg == 1
  cellsegm.show(bw3_binaryMask,5);
  CC2=imoverlay(I2,bw3_binaryMask);
  cellsegm.show(CC2,11);
  title({CurrDir+": Watershed split: "+ cellbw3stat.NumObjects + " Cells";ImgTag})
  axis off;
 % imshow(imoverlay(imsegmB,bw3,'yellow'));
  % Create yellow boundaries for the watershed overlay
  watershedBoundaryPixels = bwperim(bw3); % Boundary pixels for watershed segments
  thickWatershedBoundary = imdilate(watershedBoundaryPixels, strel('disk', 1)); % Thicken boundaries

  % Overlay the watershed boundaries onto the original image
  watershedOverlay = imoverlay(I2, thickWatershedBoundary, brightRed); % Red for watershed

  % Display the result
    imshow(watershedOverlay);
    title({CurrDir+", "+ ImgTag+": Watershed Split in Red: "+ cellbw3stat.NumObjects + " Cells"});
    axis off;
end

% Optionally save the combined overlay image
if StoreImages
    outputFolder = fullfile('C:\LCSD\Results', regexprep(CurrDir, '[\\/:*?"<>|]', ''), 'Split');
    % Clean ImgTag by replacing special characters with underscores
    cleanImgTag = regexprep(Wfile, '[\\/:*?"<>|]', '');
    % Construct the fully qualified file path
    outputFileName = fullfile(outputFolder, strcat(cleanImgTag, '.tif'));
    % Display the output filename for verification
    disp(['Saving to: ' outputFileName]);

% Create a new figure to display the image and set the title
    figureHandle = figure;
    imshow(watershedOverlay);
    title(string(CurrDir) + ", " + string(ImgTag) + ": Cell Boundaries after Split: " + num2str(cellbw3stat.NumObjects) + " Cells");

    % Save the figure with the title as a .tif image
    exportgraphics(figureHandle, outputFileName, 'Resolution', 300);

    % Close the figure to avoid clutter
    close(figureHandle);
end

BlueMask = bw3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RED channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IR = imread(theimage);
I2R = im2double(IR);
imsegmR = I2R(:,:,1); 
%figure(124); imshow(imsegmR)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Post processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D_B = bwdist(~BlueMask);
% L_B = watershed(D_B);
% L_B(~BlueMask) = 0;
CC_B = bwconncomp(BlueMask); % Finds connected components in image

% D_R = bwdist(~RedMask);
% L_R = watershed(D_R);
% L_R(~RedMask) = 0;
% CC_R = bwconncomp(RedMask);

area = []; % calculated from blue
fltarea = []; % Filtered version of area where some cells by removing very small ones
fltTrainArea = []; % Filtered version of area for training images
R_int = []; % p63
R_int_comp = []; % p63
fltR_int = []; %Filtered version of R_int
fltR_in_comp = [];
RedArea = [];
RedDensity = [];
fltRedArea = [];
fltRedDensity = [];
CellLoc = [];
fltCellLoc = []; % Filtered version of CellLoc
fltCnt = 0; %Filtered Count where very small and large cells are dropped
fltUpdtTrainCnt = 0; %Only for generating images for training
MinArea = 50; % Anything less than MinArea pixels is dropped, Default 25
MaxArea = 2000; % Anything greater than MaxArea pixels is dropped, Default 1000
fltblank = zeros(720*SF,960*SF);
fltblankMask = zeros(720*SF,960*SF);
fltTrainblank = zeros(720*SF,960*SF);
ImageSize = [720*SF 960*SF];
for k = 1:size(CC_B.PixelIdxList,2)
    cell_ind = imsegmB(CC_B.PixelIdxList{k}); %select out pixels from blue channel image

    %%%Centroid of each mask
    
    blank = zeros(720*SF,960*SF); 
    blank(CC_B.PixelIdxList{k}) = 1; % figure; imagesc(blank)
    stats = regionprops(blank,'Perimeter', 'Solidity', 'Area', 'Centroid', 'BoundingBox','Circularity', 'Solidity', 'ConvexArea', 'MajorAxisLength','MinorAxisLength','Eccentricity'); 
    CellLoc(k,1:2) = stats.Centroid;
    
    area(k) = size(CC_B.PixelIdxList{k},1); %area is # of pixel = # of indexes
    [row,col] = ind2sub(ImageSize,CC_B.PixelIdxList{k}); %Row,col index of PixelList
    CellTouchesEdge = 0;
    if ismember(1,row) || ismember(1,col) || ismember(ImageSize(1),row) || ismember(ImageSize(2),col)
        CellTouchesEdge = 1;
    end
    fltUpdt = 0;
    if area(k) > MinArea && area(k) < MaxArea
       fltUpdtTrainCnt = fltUpdtTrainCnt + 1; % This does not filter the image for Cellpose2.0
       fltTrainArea(fltUpdtTrainCnt) = area(k);
       fltTrainblank(CC_B.PixelIdxList{k}) = fltUpdtTrainCnt; % Cellpose 2.0 mask, different value for every region
  %     if (stats.Circularity > 0.7) && (stats.Solidity > 0.9) && (stats.Eccentricity <= 0.9) && (CellTouchesEdge == 0)
       if (stats.Circularity > 0.75) && (stats.Solidity > 0.95) && (stats.Eccentricity <= 0.8) && (CellTouchesEdge == 0)
        fltUpdt = 1;
        fltCnt = fltCnt + 1;
        fltarea(fltCnt) = area(k);
        fltCellLoc(fltCnt,1:2) = stats.Centroid;
        fltblank(CC_B.PixelIdxList{k}) = 255;
        fltblankMask(CC_B.PixelIdxList{k}) = fltCnt;
        CC_B.fltPixelIdxList{fltCnt}= CC_B.PixelIdxList{k};
       end
    end


    %same result could be grabbed from stats.  NIce job jun. 

    if fltUpdt
      RedCell = imsegmR(CC_B.fltPixelIdxList{fltCnt}); %%% Important! Computed (R_int_comp) red intensity is computed from red image with blue image mask
      R_int_comp(fltCnt) = mean(RedCell);
      fltR_int_comp(fltCnt) = mean(RedCell);
    end
 
    if fltUpdt
      BlueCell = imsegmB(CC_B.fltPixelIdxList{fltCnt});
      B_int_comp(fltCnt) = mean(BlueCell);
      fltB_int_comp(fltCnt) = mean(BlueCell);
    end
    %%% Note here that above red mask generation is not going to necesarily
    %%% identify the same cell
    
     if fltUpdt
      RedDensity(fltCnt) = R_int_comp(fltCnt)/size(CC_B.fltPixelIdxList{fltCnt},1);
      fltRedDensity(fltCnt) = RedDensity(fltCnt);
    end
end
if StoreForTraining
%    cellsegm.show(fltTrainblank,100)
    %saveas(gcf, WriteDir+string('\')+string('\Training\')+string(Wfile), 'tif');
    %imwrite(mat2gray(fltTrainblank), WriteDir+string('\')+string('\Training\')+string(Wfile)+string('.tif'));
    save(WriteDir+"\"+"\Training\"+string(Wfile)+".mat",'fltTrainblank')
end
if StoreForTraining1
%    cellsegm.show(fltTrainblank,100)
    %saveas(gcf, WriteDir+string('\')+string('\Training\')+string(Wfile), 'tif');
    %imwrite(mat2gray(fltTrainblank), WriteDir+string('\')+string('\Training\')+string(Wfile)+string('.tif'));
    %save(WriteDir+string('\')+string('\Filtered\')+string(Wfile)+string('.mat'),'fltblankMask')
end
if DisplayImg == 1
%  cellsegm.show(CC3,12);
%  title({CurrDir+": Shape filtered: "+ fltCnt + " Cells";ImgTag})
  cellsegm.show(fltblank,4);
  title({CurrDir+": Split cells:  " + fltCnt + " cells (Shape filter)";ImgTag})
  axis off;

  % Create red boundaries for the overlay
  filteredBoundaryPixels = bwperim(fltblank); % Boundary pixels for watershed segments
  thickFilteredBoundary = imdilate(filteredBoundaryPixels, strel('disk', 1)); % Thicken boundaries

  figureHandle = figure;
  imshow(imoverlay(I2,thickFilteredBoundary,brightRed));
  title({CurrDir+", "+ImgTag+": Shape filtered: "+ fltCnt + " Cells"});
  if StoreImages
    outputFolder = fullfile('C:\LCSD\Results', regexprep(CurrDir, '[\\/:*?"<>|]', ''), 'Filtered');
    % Clean ImgTag by replacing special characters with underscores
    cleanImgTag = regexprep(Wfile, '[\\/:*?"<>|]', '');
    % Construct the fully qualified file path
    outputFileName = fullfile(outputFolder, strcat(cleanImgTag, '.tif'));
    % Display the output filename for verification
    disp(['Saving to: ' outputFileName]);
       % Save the figure with the title as a .tif image
    exportgraphics(figureHandle, outputFileName, 'Resolution', 300);
    % Close the figure to avoid clutter
    close(figureHandle);
  end
end
cellCounts3 = fltCnt;
if StoreImageSidebySide
    if CurrDir == "CLIN\_0007" | CurrDir == "CLIN\_0011" | CurrDir == "CLIN\_0012" | CurrDir == "CLIN\_0013" | CurrDir == "CLIN\_0014"
      StatFile = statsSubDir+"_p63"+regexprep(ImgTag, '[\\/:*?"<>|]', '')+"Overlay";
      StatFullFile = statsDir+"\Analysis\"+statsSubDir+"\"+statsSubDir+"_p63\"+StatFile+".csv";
    else
      StatFile = statsSubDir+"_p63_"+regexprep(ImgTag, '[\\/:*?"<>|]', '')+"Overlay";
      StatFullFile = statsDir+"\Analysis\"+statsSubDir+"\"+statsSubDir+"_p63_\"+StatFile+".csv";
    end
    %T = readtable(StatFullFile);
    if exist(StatFullFile, 'file') == 2
        T = readtable(StatFullFile);
    else
        disp('File not found.');
    end
    fileID = fopen(StatFullFile,'r');
    C_text1 = textscan(fileID,'%s',3,'delimiter',',');
    C_text2 = textscan(fileID,'%s',3,'delimiter',',');
    textscan(fileID,'%s',1,'delimiter',',');
    textscan(fileID,'%s',9,'delimiter',',');
    C_text3 = textscan(fileID,'%s',9,'delimiter',',');
    C_text4 = textscan(fileID,'%s',9,'delimiter',',');
    C_text5 = textscan(fileID,'%s',9,'delimiter',',');
    C_text6 = textscan(fileID,'%s',9,'delimiter',',');
    C_text7 = textscan(fileID,'%s',9,'delimiter',',');
    fclose(fileID);
    CellArea = T(8:end,4);
    %imshowpair(imrotate(imsegmB,90), imrotate(rgb2gray(CC3),90), 'montage');
    %set(gca,'view',[90 90]);
    %subplot(1,2,1), imshow(imsegmB);
    %subplot(1,2,2), imshow(rgb2gray(CC3));
end
if StoreImageSidebySide1
    t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact',TileSpacing='tight' );
    nexttile
    imshow(I2);
    title(string(CurrDir) + ", " + string(ImgTag)+": Original");
    nexttile
   % imshow(rgb2gray(CC3));
    imshow(imoverlay(I2,thickFilteredBoundary,brightRed));
    title(string('Filtered, ')+ fltCnt + " Cells")
   %  saveas(gcf, WriteDir+string('\')+string('\Orig+Filtered\')+string(Wfile), 'pdf');
    outputFolder = fullfile(WriteDir, "Orig+Filtered");
    outputFileName = string(Wfile) + ".pdf";

% Ensure the directory exists before saving
    if ~exist(outputFolder, 'dir')
       mkdir(outputFolder);
    end

% Save the figure to the specified path
    saveas(gcf, fullfile(outputFolder, outputFileName));
   %   saveas(gcf, WriteDir+string('\')+string('\Orig+Filtered\')+string(Wfile), 'pdf');
 %    [counts1, binCenters1] = hist(fltarea, 25);
 %    [counts2, binCenters2] = hist(CellArea.Area, 25);
 %    figure;
 %    plot(binCenters1, counts1, 'r-');
 %    hold on;
 %    plot(binCenters2, counts2, 'g-');
 %    legend({"CellPose", "Software"});
 %    title(ImgTag+string(': Histogram'))
 %    grid on;
 % %   saveas(gcf, WriteDir+"\"+string('\Orig+Filtered\')+string(Wfile)+string('_hist'), 'pdf');
 %    saveas(gcf, WriteDir+"\"+string('\Orig+Filtered\')+string(Wfile)+"_hist", 'fig');
  %  kl = KLDiv(counts1/sum(counts1),counts2/sum(counts2)); % OneDimensional kullback liebler divergence
end

R = R_int_comp;% = cell2mat(data(12:end, 8));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Segmentation analysis - how did it do?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cd('C:\LCSD\P63 dataset 111821\P63 dataset 111821\p63 CLIN_0001\Analysis D1\D1 analysis\cLSC');
%% BZXimage = ['p63' theimage];

%figure(199)
% BZX = imread(BZXimage);
%BZX2 = im2double(BZX);

if DisplayImg == 2
  figure(5);
  hold on;
  %cellsegm.show(I,5);
  imshow(I);
  hold on;
  plot(fltCellLoc(:,1),fltCellLoc(:,2),'o', 'MarkerSize', 20, 'color', [0.8 0 0 ]); 
  title({CurrDir+": Segmented split  " + fltCnt + " cells overlaid on original image";ImgTag})
  CC=imoverlay(imsegmB,fltblank);
  figure;
  hold on;
  cellsegm.show(CC,6);
  title({CurrDir+": Segmented split  " + fltCnt + " cells overlaid on original image";ImgTag})
  axis off;
  fltcellnumber = num2cell(1:size(fltCellLoc,1));
  text(fltCellLoc(:,1),fltCellLoc(:,2),fltcellnumber,'VerticalAlignment','bottom','HorizontalAlignment','right','Color','y')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Decision time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 DensityTH = 0.45e-3;
% 
 meanTH = 85.33/255; %Threshold

% if DisplayImg == 1
%   figure(10);
%   imshow(I);
%   hold on;
% 
%   fltp63bright = 0;
%   for m = 1:fltCnt
%       if fltR_int_comp(m) > meanTH
%          plot(fltCellLoc(m,1),fltCellLoc(m,2),'o', 'MarkerSize', 20, 'color', [0 0.75 0]); 
%          text(fltCellLoc(m,1),fltCellLoc(m,2),num2str(floor(fltR_int_comp(m)*255)),'VerticalAlignment','middle','HorizontalAlignment','center','Color','g','FontSize',8)
%          fltp63bright = fltp63bright + 1;
%       else
%          plot(fltCellLoc(m,1),fltCellLoc(m,2),'o', 'MarkerSize', 20, 'color', [0.75 0 0]); 
%          text(fltCellLoc(m,1),fltCellLoc(m,2),num2str(floor(fltR_int_comp(m)*255)),'VerticalAlignment','middle','HorizontalAlignment','center','Color','r','FontSize',8)
%       end
%   end
%   title({CurrDir+": " + fltp63bright + " cells with mean > " + meanTH*255 +  "  Img:" + ImgTag + "  GreenCircle > Thresh, RedCircle < Thresh"})
% end

close all

%%% save?
save([theimage(1:end-4) '_CSoutput'])

STAT{imInx, (K-1)*4+1} = imageNames;
STAT{imInx, (K-1)*4+1 + 1} = cellCounts1;
STAT{imInx, (K-1)*4+1+ 2} = cellCounts2;
STAT{imInx, (K-1)*4+1+ 3} = cellCounts3

WriteDirTable = 'C:\LCSD\P63 dataset 111821\P63 dataset 111821';
StatTable = cell2table(STAT);

% Optional: Add headers if you have specific names for each column
% headers = ["Column1", "Column2", ..., "Column104"];
% StatTable.Properties.VariableNames = headers;

% Define file paths for saving
outputExcelFile = fullfile(WriteDirTable, 'StatData.xlsx');
outputCsvFile = fullfile(WriteDirTable, 'StatData.csv');

% Write to Excel and CSV files
writetable(StatTable, outputExcelFile, 'WriteVariableNames', false); % Excel file
writetable(StatTable, outputCsvFile, 'WriteVariableNames', false);   % CSV file

clearvars -except STAT K AnalysisList CurrDir ImageList DisplayImg StoreImages StoreImageSidebySide StoreImageSidebySide1 StoreImagesTif StoreForTraining StoreForTraining1 WriteDir statsSubDir statsDir;
end
end

function Img = checknghbr(Img)
  [maxr, maxc] = size(Img);
  for r = 1:maxr
     for c = 1:maxc
        if Img(r,c) ~= 0
          if r-1 >= 1 && c-1 >= 1  % Top Left 
             if Img(r-1,c-1) ~= Img(r,c) && Img(r-1,c-1) ~= 0
                Img(r,c) = 0;
             end
          end
          if r-1 >= 1   % Top middle
             if Img(r-1,c) ~= Img(r,c) && Img(r-1,c) ~= 0
                Img(r,c) = 0;
             end
          end
          if r-1 >= 1 && c+1 <= maxc  % Top right
             if Img(r-1,c+1) ~= Img(r,c) && Img(r-1,c+1) ~= 0
                Img(r,c) = 0;
             end
          end
          if  c-1 >= 1  % Left same row
             if Img(r,c-1) ~= Img(r,c) && Img(r,c-1) ~= 0
                Img(r,c) = 0;
             end
          end
          if c+1 <= maxc  % Right same row
             if Img(r,c+1) ~= Img(r,c) && Img(r,c+1) ~= 0
                Img(r,c) = 0;
             end
          end
          % if r+1 <= maxr && c-1 >= 1  % Bottom left
          %    if Img(r+1,c-1) ~= Img(r,c) && Img(r+1,c-1) ~= 0
          %       Img(r,c) = 0;
          %    end
          % end
          if r+1 <= maxr   % Bottom middle
             if Img(r+1,c) ~= Img(r,c) && Img(r+1,c) ~= 0
                Img(r,c) = 0;
             end
          end
          % if r+1 <= maxr && c+1 <= maxc  % Bottom right
          %    if Img(r+1,c+1) ~= Img(r,c) && Img(r+1,c+1) ~= 0
          %       Img(r,c) = 0;
          %    end
          % end
        end
     end
  end
end
function dist=KLDiv(P,Q)
%  dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
%  distributions
%  P and Q  are automatically normalised to have the sum of one on rows
% have the length of one at each 
% P =  n x nbins
% Q =  1 x nbins or n x nbins(one to one)
% dist = n x 1
if size(P,2)~=size(Q,2)
    error('the number of columns in P and Q should be the same');
end
if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
   error('the inputs contain non-finite values!') 
end
% normalizing the P and Q
if size(Q,1)==1
    Q = Q ./sum(Q);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./repmat(Q,[size(P,1) 1]));
    temp(isnan(temp))=0;% resolving the case when P(i)==0
    dist = sum(temp,2);
    
    
elseif size(Q,1)==size(P,1)
    
    Q = Q ./repmat(sum(Q,2),[1 size(Q,2)]);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./Q);
    temp(isnan(temp))=0; % resolving the case when P(i)==0
    dist = sum(temp,2);
end
end