clear; close all; clc;

%% Parameters and Initialization
StoreTiledImages = 1;
SaturationVal = 0.7;
remove_last_point = 1;
    % Applies the selected thresholding scheme to the input data.
    % Valid choices for 'scheme' are:
    % 'Otsu', 'Li', 'Kapur', 'IsoData', 'Entropy', 'Renyi', 'Kittler', 'Shanbhag'
schemes = {'Otsu', 'Li', 'Kapur', 'IsoData'};

% Directory and file setup
cd("C:\LCSD\P63 dataset 111821");
ExcelFname = 'AnalysisBookForSegmentationCompare.xlsx';
AnalysisList = readtable(ExcelFname, 'Format', 'auto');
DataIdx1 = 1;

% Histogram settings
nbins = 255;
Histedges = linspace(0, 1, nbins);

% Plot settings
CircleSegmentatedCells = 0;

% Predefine variables to avoid resizing inside loops.
maxFiles = 100; % Adjust based on expected maximum number of files.
AveT = {};
% Initialize storage for thresholds
thresholds = struct();
%% Loop Through Each Analysis Folder
for K = 1:size(AnalysisList.Folder, 1)
    % Extract experiment parameters
    DengThreshold = AnalysisList.DengThreshold(K);
    CurrentFolder = AnalysisList.Folder{K};
    disp(CurrentFolder)
    ExpTitle = AnalysisList.Experiment{K};
    ThreshD = AnalysisList.ThreshD{K};
    ThreshR = AnalysisList.ThreshR{K};
    Pdir = AnalysisList.PDir{K};
    cd(CurrentFolder);
    GT1 = struct(); % Empty structure to store threshold types
    GT2 = struct();
    GT3 = struct();
    GT4 = struct();
    % List and sort data files
    DataList = ls('*Overlay*.mat');
    DataList = natsortfiles(DataList);
    if Pdir == "P15" && ThreshD == "D1"
        RnoPrev = char('3'); % Special case of P0015, D1 where R's begin with 3 (there is no R1 and R2)
    else
        RnoPrev = char('1');
    end
    % Initialize variables
    DataIdxLast = 1;
    RfileNo = 0;
    DataIdx1Memory = DataIdx1;

    %% Process Each File
    for DataIdx = 1:size(DataList, 1)
        % Load data
        fname = strtrim(DataList(DataIdx, :));
        load(fname);
        DataIdx1 = DataIdx1Memory;
        STAT = cell(size(DataList,1), 12); 
        wfname = strrep(fname, '_Overlay_CSoutput.mat', '');
        Rno = fname(strfind(fname,'R') + 1); % find the number next to R in filename

        % Update average flag and R-number check
        UpdateAverage = (Rno ~= RnoPrev || DataIdx == size(DataList, 1));
        if Rno ~= RnoPrev, RfileNo = 1; else, RfileNo = RfileNo + 1; end
        RnoPrev = Rno;

        %% Red Cell Intensity Calculations
        imsegmR_filtered = imsegmR .* fltblank / 255.0;
        imsegmR_limited = min(imsegmR_filtered, SaturationVal);
        [RedMean, RedMean_70] = calculate_cell_means(CC_B, imsegmR_filtered, imsegmR_limited, remove_last_point, SaturationVal);

        %% Segmentation and Thresholding
        addpath 'C:\LCSD\P63 dataset 111821';
         % Perform thresholding
        [segmentedCells, cellStats, thresholds] = optimized_segmentation_and_thresholding( ...
            imsegmR, fltblank, SaturationVal, remove_last_point, Histedges, nbins, schemes);


        % Store all thresholds dynamically
        [GT1,GT2,GT3,GT4] = storeAllThresholdResults(thresholds, DataIdx, GT1,GT2,GT3,GT4);
        %% Update Averages for Current Experiment
        if UpdateAverage
          [AveT1, AveT2, AveT3, AveT4, AveT5, AveT, DataIdxLast, DataIdx1] = updateAverageValues(DataIdx, DataList, DataIdxLast, ...
                                                                                  DataIdx1, DengThreshold, GT1, GT2, GT3, GT4, ...
                                                                                  ThreshD, Rno, ThreshR, Pdir, AveT);
        end
    end

    close all
    clearvars -except K STAT AveT schemes Histedges nbins SaturationVal ExcelFname DengThreshold remove_last_point DataIdx1 ForPlots StoreTiledImages
    AnalysisList = readtable(ExcelFname,'Format','auto');
end
plotThresholdComparisons(schemes, AveT);


%% Helper Functions
function [RedMean, RedMean_70] = calculate_cell_means(CC_B, imsegmR, imsegmR_limited, remove_last_point, SaturationVal)
    fltCnt = CC_B.NumObjects;
    RedMean = zeros(1, fltCnt);
    RedMean_70 = zeros(1, fltCnt);

    for zz = 1:fltCnt
        pixelIdxs = CC_B.PixelIdxList{zz};
        pixelValues = imsegmR(pixelIdxs);
        pixelValues_70 = imsegmR(pixelIdxs);

        if remove_last_point
            pixelValues = pixelValues(pixelValues ~= 1.0);
            pixelValues_70 = pixelValues_70(pixelValues_70 ~= 1.0);
        end
        pixelValues_70(pixelValues_70>SaturationVal) = SaturationVal;

        if ~isempty(pixelValues), RedMean(zz) = mean(pixelValues); end
        if ~isempty(pixelValues_70), RedMean_70(zz) = mean(pixelValues_70); end
        
    end
    RedMean = RedMean(RedMean > 0);
    RedMean_70 = RedMean_70(RedMean_70 > 0);
end

function [segmentedCells, cellStats, thresholds] = optimized_segmentation_and_thresholding( ...
    imsegmR, fltblank, SaturationVal, remove_last_point, Histedges, nbins, schemes)

    % Apply segmentation and thresholding
    imsegmR1 = imsegmR.*fltblank/255.0; % Important, filter out the Red cells with the Blue mask

    % Step 1: Attenuate the red channel image
    imsegmR_attenuated = imsegmR1;
    imsegmR_attenuated(imsegmR_attenuated > SaturationVal) = SaturationVal;

    % Step 2: Segment individual cells (Connected Components Analysis)
    CC_B = bwconncomp(fltblank);
    fltCnt = CC_B.NumObjects;

    % Initialize statistics arrays
    RedMean = zeros(1, fltCnt);
    RedMean_70 = zeros(1, fltCnt);

    % Compute mean intensity per segment
    for zz = 1:fltCnt
        pixelIdxs = CC_B.PixelIdxList{zz};
        pixelValues = imsegmR(pixelIdxs);

        % Exclude pixel values at limits if necessary
        if remove_last_point
            pixelValues = pixelValues(pixelValues ~= 1.0);
        end
        pixelValues_70 = min(pixelValues, SaturationVal);

        % Compute mean intensity values
        if ~isempty(pixelValues)
            RedMean(zz) = mean(pixelValues);
            CellArea(zz) = length(pixelValues);
        end
        if ~isempty(pixelValues_70)
            RedMean_70(zz) = mean(pixelValues_70);
        end
    end

    % Histogram-based thresholding
    [counts1, edges1] = imhist(imsegmR1, nbins);
    counts1(1) = 0; % Exclude the background pixel value
    if remove_last_point
        counts1 = counts1(1:end-1);
        edges1 = edges1(1:end);
    end
    
    if remove_last_point
        imsegmR_temp = imsegmR1(imsegmR1 ~= 1.0);
        imsegmR_temp = min(imsegmR_temp, SaturationVal);
    else 
        imsegmR_temp = min(imsegmR1, SaturationVal);
    end
    [counts2, edges2] = imhist(imsegmR_temp, nbins);
    counts2(1) = 0; % Exclude the background pixel value

    histo_meanR = histcounts(RedMean, Histedges);
    histo_meanR_70 = histcounts(RedMean_70, Histedges);

    % Now computing thresholds
    %global_threshold = applyThreshold(counts1, edges1, scheme);
    %attenuated_threshold = applyThreshold(counts2, edges2, scheme);

    % Per-cell histogram thresholding
    % mean_threshold = applyThreshold(histo_meanR, Histedges, scheme);
    %mean_70_threshold = applyThreshold(histo_meanR_70, Histedges, scheme);
    % CLEANUP
    % Compute thresholds for each scheme
    for i = 1:length(schemes)
      scheme = schemes{i};
      switch scheme
      case 'Otsu'
        thresholds.regular_threshold.(scheme) = otsuthresh(counts1);
        thresholds.attn_threshold.(scheme) = otsuthresh(counts2);
        thresholds.mean_threshold.(scheme) = otsuthresh(histo_meanR);
        thresholds.meanattn_threshold.(scheme) = otsuthresh(histo_meanR_70);
      case 'Li'
        thresholds.regular_threshold.(scheme) = li_threshold(counts1, edges1);
        thresholds.attn_threshold.(scheme) = li_threshold(counts2, edges2);
        thresholds.mean_threshold.(scheme) = li_threshold(histo_meanR, Histedges);
        thresholds.meanattn_threshold.(scheme) = li_threshold(histo_meanR_70, Histedges);
      case 'Kapur'
        thresholds.regular_threshold.(scheme) = kapur_threshold(counts1, edges1);
        thresholds.attn_threshold.(scheme) = kapur_threshold(counts2, edges2);
        thresholds.mean_threshold.(scheme) = kapur_threshold(histo_meanR, Histedges);
        thresholds.meanattn_threshold.(scheme) = kapur_threshold(histo_meanR_70, Histedges);
      case 'IsoData'
        thresholds.regular_threshold.(scheme) = isodata_threshold(counts1, edges1);
        thresholds.attn_threshold.(scheme) = isodata_threshold(counts2, edges2);
        thresholds.mean_threshold.(scheme) = isodata_threshold(histo_meanR, Histedges);
        thresholds.meanattn_threshold.(scheme) = isodata_threshold(histo_meanR_70, Histedges);
      case 'Entropy'
        thresholds.regular_threshold.(scheme) = entropyThresholding(counts1, edges1); % Example alpha = 2
        thresholds.attn_threshold.(scheme) = entropyThresholding(counts2, edges2);
        thresholds.mean_threshold.(scheme) = entropyThresholding(histo_meanR, Histedges);
        thresholds.meanattn_threshold.(scheme) = entropyThresholding(histo_meanR_70, Histedges);
      otherwise
        error('Unknown ThresholdScheme: %s', scheme);
      end
    end

    % Populate outputs
    thresholds = struct('regular_threshold', thresholds.regular_threshold, ...
                        'attn_threshold', thresholds.attn_threshold, ...
                        'mean_threshold', thresholds.mean_threshold, ...
                        'mean_attn_threshold', thresholds.meanattn_threshold);

    segmentedCells = CC_B; % Segmented cell data
    cellStats = struct('RedMean', RedMean, 'RedMean_70', RedMean_70);
end


function [AveT1, AveT2, AveT3, AveT4, AveT5, AveT, DataIdxLast, DataIdx1] =  updateAverageValues(DataIdx, DataList, DataIdxLast, ...
                                                                                   DataIdx1, DengThreshold,T1, T2, T3, T4,  ...
                                                                                   ThreshD, Rno, ThreshR, Pdir, AveT)
% Define threshold field names (ensure this is available beforehand)
thresholdFields = fieldnames(T1); % Extract all threshold types from T dynamically

if DataIdx == size(DataList, 1)
    % Calculate averages and update AveT for the last DataIdx
    [AveT1, AveT2, AveT3, AveT4, AveT5] = calculateThresholds(DataIdx, DataIdxLast, T1, T2, T3, T4, thresholdFields, DengThreshold);
    AveT = transferDataToAveT(DataIdx, DataIdx1, AveT1, AveT2, AveT3, AveT4, AveT5, ThreshD, strcat('R', char(Rno)), ThreshR, Pdir, AveT);
    DataIdxLast = DataIdx + 1;
    DataIdx1 = DataIdx1 + 1;
else
    % Calculate averages and update AveT for intermediate DataIdx
    [AveT1, AveT2, AveT3, AveT4, AveT5] = calculateThresholds(DataIdx - 1, DataIdxLast, T1, T2, T3, T4, thresholdFields, DengThreshold);
    AveT = transferDataToAveT(DataIdx - 1, DataIdx1, AveT1, AveT2, AveT3, AveT4, AveT5, ThreshD, strcat('R', char(Rno)), ThreshR, Pdir, AveT);
    DataIdxLast = DataIdx;
    DataIdx1 = DataIdx1 + 1;
end

end

function threshold = li_threshold(counts, edges)
    % Compute histogram bin centers
    binCenters = (edges(1:end-1) + edges(2:end)) / 2;

    % Normalize histogram with small epsilon
    epsilon = 1e-12;
    counts = counts + epsilon;
    counts = counts / sum(counts);

    % Initial guess for threshold
    T_prev = mean(binCenters);

    % Iterative optimization
    while true
        % Divide bins into two groups based on the current threshold
        g1 = binCenters(binCenters <= T_prev);
        g2 = binCenters(binCenters > T_prev);

        idx1 = binCenters <= T_prev;
        idx2 = binCenters > T_prev;

        p1 = sum(counts(idx1));
        p2 = sum(counts(idx2));

        % Compute means for both groups
        m1 = sum(g1 .* counts(idx1)) / p1;
        m2 = sum(g2 .* counts(idx2)) / p2;

        % Update threshold
        T_new = (m1 + m2) / 2;

        % Convergence check
        if abs(T_prev - T_new) < 1e-6
            break;
        end
        T_prev = T_new;
    end

    % Final threshold
    threshold = T_prev;
end

function th = kapur_threshold(counts, edges)
    % Normalize histogram
    prob = counts / sum(counts);
    GradeI = length(prob);  % Dynamically determine histogram length
    psai = zeros(GradeI, 1);    
    
    prob_t = 0;
    entropy_t = 0;              
    
    % Total entropy
    ind = prob > 0;
    entropy_L = -sum(prob(ind) .* log(prob(ind) + eps));
    
    for i = 1:GradeI
        prob_t = prob_t + prob(i);
        
        if prob(i) > 0 && prob_t < 1
            entropy_t = entropy_t - prob(i) * log(prob(i) + eps);
            psai(i) = log(prob_t * (1 - prob_t) + eps) + ...
                      entropy_t / (prob_t + eps) + ...
                      (entropy_L - entropy_t) / (1 - prob_t + eps);
        elseif prob(i) == 0 && i > 1
            psai(i) = psai(i - 1);  % Copy previous value
        end
    end
    
    [~, ind] = max(psai);
    th = edges(ind);  % Return threshold value
end

function threshold = isodata_threshold(counts, edges)
    % ISODATA_THRESHOLD Computes the ISODATA threshold using histogram data.
    % INPUT:
    %   counts: Histogram bin counts (array)
    %   edges: Histogram bin edges (array)
    % OUTPUT:
    %   threshold: Computed threshold value (scalar)

    % Step 1: Calculate bin centers from edges
    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

    % Align `counts` with `bin_centers` by removing the last bin (if necessary)
    if length(counts) == length(bin_centers) + 1
        counts = counts(1:end-1);
    end

    % Ensure `counts` is a row vector (to match `bin_centers`)
    counts = counts(:)'; % Convert to row vector

    % Step 2: Expand histogram counts to represent data points
    try
        data = repelem(bin_centers, counts);
    catch
        error('Error in repelem: counts must contain non-negative integers.');
    end

    % Remove zeros (background values), if necessary
    data = data(data > 0);

    % Step 3: Initialize the threshold (use the mean of the data)
    t_prev = mean(data);
    tol = 0.01;          % Convergence tolerance
    max_iter = 100;      % Maximum iterations

    % Step 4: Iterative ISODATA thresholding
    for iter = 1:max_iter
        % Separate data into two groups based on the current threshold
        group1 = data(data <= t_prev);
        group2 = data(data > t_prev);

        % If one of the groups is empty, stop (threshold cannot improve)
        if isempty(group1) || isempty(group2)
            break;
        end

        % Compute the means of the two groups
        mean1 = mean(group1);
        mean2 = mean(group2);

        % Update the threshold as the mean of the two group means
        t_new = (mean1 + mean2) / 2;

        % Check for convergence
        if abs(t_new - t_prev) < tol
            break;
        end

        % Update threshold for the next iteration
        t_prev = t_new;
    end

    % Step 5: Output the final threshold
    threshold = t_prev;
end

function [AveT1, AveT2, AveT3, AveT4, AveT5] = calculateThresholds(CurrentIdx, LastIdx, T1, T2, T3, T4, thresholdFields, DengThreshold)
    % calculateThresholds Computes average values for thresholds
    %
    % Inputs:
    %   CurrentIdx       - Current data index
    %   LastIdx          - Last processed data index
    %   T                - Structure with all thresholds
    %   thresholdFields  - List of threshold field names
    %   DengThreshold    - Deng threshold value
    %
    % Outputs:
    %   T1, T2, T3, T4, T5 - Ordered average values of thresholds

    numFields = min(length(thresholdFields), 6); % Support up to 6 dynamic thresholds
    thresholdAverages = nan(1, numFields);

    % Loop through the threshold fields and compute averages
    for i = 1:numFields
        fieldName = thresholdFields{i};
        thresholdValues1 = T1.(fieldName)(LastIdx:CurrentIdx);
        thresholdAverages1(i) = nanmean(thresholdValues1); % Use nanmean to handle NaNs
        thresholdValues2 = T2.(fieldName)(LastIdx:CurrentIdx);
        thresholdAverages2(i) = nanmean(thresholdValues2); % Use nanmean to handle NaNs
        thresholdValues3 = T3.(fieldName)(LastIdx:CurrentIdx);
        thresholdAverages3(i) = nanmean(thresholdValues3); % Use nanmean to handle NaNs
        thresholdValues4 = T4.(fieldName)(LastIdx:CurrentIdx);
        thresholdAverages4(i) = nanmean(thresholdValues4); % Use nanmean to handle NaNs
    
        
        % Assign to T1, T2, ..., T5 in order
        AveT1.(fieldName) = thresholdAverages1(i); % First threshold average
        AveT2.(fieldName) = thresholdAverages2(i); % Second threshold average
        AveT3.(fieldName) = thresholdAverages3(i); % Third threshold average
        AveT4.(fieldName) = thresholdAverages4(i); % Fourth threshold average
        AveT5.(fieldName) = DengThreshold;        % DengThreshold assigned to T5
    end
end
function AveT = transferDataToAveT(currentIdx, DataIdx1, AveT1, AveT2, AveT3, AveT4, AveT5, ThreshD, RnoStr, ThreshR, Pdir, AveT)
    % Transfer calculated thresholds and metadata to AveT

    % Combine AveT1 to AveT5 into a cell array for easier processing
    AveTStructs = {AveT1, AveT2, AveT3, AveT4};
    
    % Get the field names (assuming all structures share the same fields)
    fieldNames = fieldnames(AveT1);

    % Initialize a temporary cell array for threshold values
    tempThresholds = cell(1, numel(AveTStructs));
    
    % Iterate over each structure (AveT1 to AveT5)
    for sIdx = 1:numel(AveTStructs)
        currentStruct = AveTStructs{sIdx};
        
        % Initialize a temporary array to store rounded values
        tempValues = zeros(1, numel(fieldNames));
        
        % Loop through each field
        for fIdx = 1:numel(fieldNames)
            fieldName = fieldNames{fIdx};
            tempValues(fIdx) = round(currentStruct.(fieldName) * 255, 3);
        end
        
        % Store the rounded values in the tempThresholds array
        tempThresholds{sIdx} = tempValues;
    end

    % Assign the calculated thresholds to AveT
    for sIdx = 1:numel(tempThresholds)
        AveT{sIdx, DataIdx1} = tempThresholds{sIdx};
    end

    % Handle additional metadata
    AveT{5, DataIdx1} = AveT5.('Otsu');  % DengThreshold assigned to T5
    AveT{6, DataIdx1} = ThreshD;         % Deng Threshold Dvalue
    AveT{7, DataIdx1} = ThreshR;         % Region Threshold
    AveT{8, DataIdx1} = Pdir;            % Parameter Direction
end

function [T1,T2,T3,T4] = storeAllThresholdResults(thresholds, DataIdx, T1,T2,T3,T4)
% storeAllThresholdResults Captures all threshold types dynamically
%   
% Inputs:
%   thresholds   - Structure with multiple threshold types (each containing subfields)
%   DataIdx      - Current index for data storage
%   T1-T4        - Structure to store threshold values
%
% Outputs:
%   T            - Updated structure with stored threshold values

    thresholdNames = fieldnames(thresholds); % Top-level threshold names
    fieldName1 = thresholdNames{1}; % Current threshold type (e.g., 'regular_threshold')
    fieldName2 = thresholdNames{2};
    fieldName3 = thresholdNames{3};
    fieldName4 = thresholdNames{4};
    subFieldNames = fieldnames(thresholds.(fieldName1)); % Subfields like 'Otsu', 'Li', etc.
    %T1.(subFieldNames) = NaN(1, 100); 
    for j = 1:length(subFieldNames)
                subFieldName = subFieldNames{j}; % Current subfield (e.g., 'Otsu')
                T1.(subFieldName)(DataIdx) = thresholds.(fieldName1).(subFieldName);
                T2.(subFieldName)(DataIdx) = thresholds.(fieldName2).(subFieldName);
                T3.(subFieldName)(DataIdx) = thresholds.(fieldName3).(subFieldName);
                T4.(subFieldName)(DataIdx) = thresholds.(fieldName4).(subFieldName);
    end
end
function optimalThreshold = entropyThresholding(counts, edges)
    % ENTROPYTHRESHOLDING: Computes optimal threshold using Shannon Entropy
    %
    % counts: Histogram counts (output of imhist)
    % edges: Bin edges corresponding to the histogram (output of imhist)
    % Returns:
    % - optimalThreshold: Calculated threshold value (bin edge)

    % Normalize the histogram counts to get a probability distribution
    totalPixels = sum(counts);  % Total number of pixels
    prob = counts / totalPixels;  % Probability distribution

    % Initialize variables to store maximum entropy
    maxEntropy = -Inf;
    optimalThreshold = 0;

    % Add a small constant to avoid log(0) and ensure positive probabilities
    eps = 1e-10;

    % Loop through all possible thresholds (1 to length-1 to split into two classes)
    for t = 1:length(prob)-1
        % Split the histogram into two classes
        prob1 = prob(1:t);
        prob2 = prob(t+1:end);

        % Calculate the entropy for both classes
        entropy1 = -sum(prob1 .* log2(prob1 + eps)); % Add eps to avoid log(0)
        entropy2 = -sum(prob2 .* log2(prob2 + eps)); % Add eps to avoid log(0)

        % Compute the total entropy
        totalEntropy = entropy1 + entropy2;

        % Check if the current entropy is the maximum found
        if totalEntropy > maxEntropy
            maxEntropy = totalEntropy;
            optimalThreshold = t;  % Save the index of the threshold
        end
    end

    % The threshold value is the bin edge corresponding to the index of the optimal threshold
    optimalThreshold = edges(optimalThreshold);
end

function plotThresholdComparisons(schemes, AveT)
    % Number of rows (thresholds) and columns (schemes)
    numRows = 5; % The first 4 rows contain threshold values
    numCols = size(AveT, 2); % Number of schemes is determined by columns in AveT
    
    % Check if AveT has enough rows for thresholds and DengThresholds
    if size(AveT, 1) < 5
        error('AveT does not have enough rows to access row 5 for DengThresholds.');
    end
    
    % Extract threshold values and DengThresholds from AveT
    thresholdRows = AveT(1:numRows, :); % Rows 1 to 4 contain threshold values
    dengThresholds = cell2mat(AveT(5, :)); % Row 5 contains DengThresholds for each scheme
    
    % Iterate through the rows (1 to 4) for threshold comparison
    numSchemes = length(thresholdRows{1,1});
    Row1 = thresholdRows(1, :);
    Row2 = thresholdRows(2, :);
    Row3 = thresholdRows(3, :);
    Row4 = thresholdRows(4, :);
    for schemeIdx = 1:numSchemes
        for colIdx = 1:numCols
            ThrVal1(schemeIdx, colIdx) = Row1{1, colIdx}(schemeIdx);
            ThrVal2(schemeIdx, colIdx) = Row2{1, colIdx}(schemeIdx);
            ThrVal3(schemeIdx, colIdx) = Row3{1, colIdx}(schemeIdx);
            ThrVal4(schemeIdx, colIdx) = Row4{1, colIdx}(schemeIdx);
        end
    end
    schemes = cat(2, schemes, {'Deng'});

    % Concatenate rows 8, 7, and 6 from AveT for x-axis labels
    row6 = string(AveT(6, :)); % Convert row 6 to string
    row7 = string(AveT(7, :)); % Convert row 7 to string
    row8 = string(AveT(8, :)); % Convert row 8 to string

    % Combine them into a single label array without a separator
    xLabels = strcat(row8, row7, row6);

    for plotIdx = 1:4
        figure;

        thrVal = eval(sprintf('ThrVal%d', plotIdx)); % Get ThrVal1, ThrVal2, etc.

        % Initialize array to store correlation values
        correlations = zeros(numSchemes, 1);

        % Plot each threshold scheme and compute the correlation with DengThresholds
        for z = 1:numSchemes
            plot(1:numCols, thrVal(z, :), '-o', 'DisplayName', schemes{z});
            hold on;
            
            % Compute the correlation between the current scheme's threshold values and DengThresholds
            correlations(z) = corr(thrVal(z, :)', dengThresholds(:));  % Transpose threshold values for correct correlation
        end

        % Plot Deng thresholds
        plot(1:length(dengThresholds), dengThresholds, '--o', ...
            'DisplayName', 'Deng Threshold', 'LineWidth', 1.5);
        hold off;

        % Customize the plot
        xlabel('Image (P,D,R) Number');
        ylabel('Threshold Values');
        set(gca, 'XTick', 1:numCols, 'XTickLabel', xLabels, 'XTickLabelRotation', 90); % Vertical labels
% Add the legend and show it
lgd = legend('Location', 'northeastoutside'); % Place legend outside
legend show;

% Get the position of the legend
legendPosition = lgd.Position; % Position of the legend in normalized coordinates

% Adjust the starting position for the text annotations below the legend
textStartX = 1.31*(legendPosition(1) + legendPosition(3)); % Align with the center of the legend
textStartY = legendPosition(2) - 0.05; % Slightly below the bottom of the legend (adjust as needed)

% Convert the plot position to normalized coordinates
axesPosition = gca().Position; % Correctly accessing the axes position

% Place the correlation text under the legend
for z = 1:numSchemes
    % Adjust vertical spacing for the annotations
    annotationY = textStartY - (z * 0.05); % Adjust vertical offset as needed
    text(textStartX, annotationY, ...
        ['Cor (' schemes{z} '): ' num2str(correlations(z), '%.2f')], ...
        'Units', 'normalized', ... % Use normalized coordinates for consistent placement
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 10, ...
        'Color', 'black');
end

        % Title customization
        switch plotIdx
            case 1
                title(['Regular Threshold (Red Cell ' ...
                    'Pixel Intensity)']);
            case 2
                title('Saturated Threshold (70% Saturation of Red Cell Pixel Intensity)');
            case 3
                title('Threshold on Mean of Red Cell Values ');
            case 4
                title('Threshold on Mean of 70% Red Cell Values');
        end

        grid on;
    end
end
