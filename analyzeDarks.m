function darkAnalysis = analyzeDarks(folder, zeroResults, gain, plotFlag)

if nargin < 3|| isempty(gain)
    gain = 1; % Assume 1 e-/ADU if not provided
end
if nargin < 4
    plotFlag = 0;
end

exposureTime = 3600;
transientThreshold = 1000;

darkStack = scripts.getImageStack(folder);
[rows, cols, numFrames] = size(darkStack);

darksCorrected = double(darkStack) - zeroResults.masterZero;

masterDark = median(darksCorrected, 3);
masterDark(masterDark < 0) = 0;

masterDark_e = masterDark * gain;
darkCurrentMap_e_s = masterDark_e / exposureTime;


% Identify Static Hot Pixels (Sensor Defects)
% Hot Pixel is a permanent sensor defect (High in ALL frames -> High Median)
meanDark = median(darkCurrentMap_e_s(:));
stdDarkSpatial = std(darkCurrentMap_e_s(:));
hotPixelThreshold = meanDark + (5 * stdDarkSpatial);
isHotPixel = darkCurrentMap_e_s > hotPixelThreshold;
numHotPixels = sum(isHotPixel(:));

kernel = ones(5, 5);
neighborCount = conv2(double(isHotPixel), kernel, 'same');

% Mask for isolated pixels: Must be Hot ITSELF (1) and Neighbor Sum must be 1
isIsolated = isHotPixel & (neighborCount == 1);
numIsolated = sum(isIsolated(:));
numClustered = numHotPixels - numIsolated;

% Transient Events (Cosmic Rays)
maxFrame = max(darksCorrected, [], 3);
isCosmicRay = (maxFrame - masterDark) > transientThreshold;
numCosmicRays = sum(isCosmicRay(:));

% 5. Transient Cleaning
darksCleaned = darksCorrected;
for i = 1:numFrames
    thisFrame = darksCorrected(:,:,i);
    diffMap = abs(thisFrame - masterDark);
    isTransient = diffMap > transientThreshold;
    thisFrame(isTransient) = masterDark(isTransient);
    darksCleaned(:,:,i) = thisFrame;
end

% 6. Noise Analysis
observedNoise = std(darksCorrected, 0, 3);
observedNoiseCleaned = std(darksCleaned, 0, 3);

% 7. Pack Results
darkAnalysis.masterDark = masterDark;
darkAnalysis.darkCurrentMap = darkCurrentMap_e_s;
darkAnalysis.observedNoiseMap = observedNoise;
darkAnalysis.observedNoiseCleanedMap = observedNoiseCleaned;
darkAnalysis.hotPixelCount = numHotPixels;
darkAnalysis.hotPixelsIsolated = numIsolated;
darkAnalysis.hotPixelsClustered = numClustered;
darkAnalysis.meanDarkCurrent = mean(darkCurrentMap_e_s(:));
darkAnalysis.meanDarkNoise = mean(observedNoise(:));
darkAnalysis.meanDarkNoiseCleaned = mean(observedNoiseCleaned(:));

% Console Output
fprintf('------------------------------------------------\n');
fprintf('DARK FRAMES ANALYSIS RESULTS (16-bit)\n');
fprintf('------------------------------------------------\n');
fprintf('Analysis Configuration:\n');
fprintf('  Frames:            %d\n', numFrames);
fprintf('  Exposure Time:     %.1f s\n', exposureTime);
fprintf('  Gain:              %.2f e-/ADU\n', gain);
fprintf('------------------------------------------------\n');
fprintf('Master Dark Current Statistics:\n');
fprintf('  Master Dark Std:   %.2f ADU (Spatial)\n', std(masterDark(:)));
fprintf('  Master Dark Mean:  %.2f ADU (Spatial)\n', mean(masterDark(:)));
fprintf('------------------------------------------------\n');
fprintf('Sensor Statistics:\n');
fprintf('  Std Dark Current:  %.7f e-/s (Uniformity)\n', stdDarkSpatial);
fprintf('  Mean Dark Current: %.7f e-/s/pixel\n', darkAnalysis.meanDarkCurrent);
fprintf('  Mean Dark Noise:   %.2f ADU (RMS)\n', darkAnalysis.meanDarkNoise);
fprintf('  Mean Dark Noise:   %.2f ADU (RMS - Cleaned)\n', darkAnalysis.meanDarkNoiseCleaned);
fprintf('------------------------------------------------\n');
fprintf('Defects & Radiation:\n');
fprintf('  Static Hot Pixels: %d (%.3f%% of sensor)\n', numHotPixels, (numHotPixels/(rows*cols))*100);
fprintf('    - Isolated:      %d (%.1f%% of hot pixels)\n', numIsolated, (numIsolated/numHotPixels)*100);
fprintf('    - Clustered:     %d (Likely sensor etching defects)\n', numClustered);
fprintf('  Gamma/Cosmic Hits: %d (Transient Events)\n', numCosmicRays);
fprintf('------------------------------------------------\n\n');
% 7. Plotting & Export
if plotFlag
    % Check/Create Export Directory
    if ~exist('img', 'dir')
        mkdir('img');
    end

    % Config
    figWidth = 1000; 
    figHeight = 800;
    mainFontSize = 18;

    % --- FIGURE 1: Master Dark ---
    f1 = figure('Name', 'Master Dark', 'Color', 'w', 'Units', 'pixels', ...
        'Position', [100 100 figWidth figHeight]);
    
    imagesc(darkAnalysis.masterDark);
    colorbar;
    colormap(gca, 'gray'); 
    % title('Master Dark (Gamma Removed)', 'FontSize', mainFontSize);
    clim([0 prctile(masterDark(:), 99)]); 
    set(gca, 'FontSize', mainFontSize);
    axis image;

    exportgraphics(f1, fullfile('img', 'dark_master.png'), 'Resolution', 300);

    % --- FIGURE 2: Dark Current Distribution (Poisson) ---
    f2 = figure('Name', 'Dark Distribution', 'Color', 'w', 'Units', 'pixels', ...
        'Position', [150 150 figWidth figHeight]);
    
    hold on;
    data_dark = darkAnalysis.darkCurrentMap(:);
    
    histogram(data_dark, 200, 'Normalization', 'pdf', 'EdgeColor', 'none');
    
    % Poisson Model
    lambda_total = meanDark * exposureTime; 
    x_rate = linspace(min(data_dark), prctile(data_dark, 99.9), 500);
    x_counts = round(x_rate * exposureTime);
    y_prob = poisspdf(x_counts, lambda_total);
    y_fit = y_prob * exposureTime;
    
    plot(x_rate, y_fit, 'r-', 'LineWidth', 3);
    xlim([0, 0.01]);
    % set(gca, 'YScale', 'log'); 
    %title(sprintf('Dark Current Distribution\n\\mu=%.7f ADU/s, \sigma=%.7f ADU/s, (Poisson Fit)', meanDark, stdDarkSpatial), ...
    %    'FontSize', mainFontSize);

    title(sprintf('\\mu=%.7f e-/s, \\sigma=%.7f e-/s, (Poisson Fit)', meanDark, stdDarkSpatial), ...
        'FontSize', mainFontSize);
    xlabel('Dark Current (e-/s)', 'FontSize', mainFontSize); 
    ylabel('Probability Density', 'FontSize', mainFontSize);
    legend({'Observed Data', 'Poisson Model'}, 'FontSize', mainFontSize, 'Location', 'best');
    set(gca, 'FontSize', mainFontSize);
    grid on;
    hold off;

    exportgraphics(f2, fullfile('img', 'dark_current_distribution.png'), 'Resolution', 300);

    % --- FIGURE 3: Hot Pixel Map (Updated with Clusters) ---
    f3 = figure('Name', 'Hot Pixel Map', 'Color', 'w', 'Units', 'pixels', ...
        'Position', [200 200 figWidth figHeight]);
    
    % Create a visualization map: 0=Normal, 1=Isolated, 2=Clustered
    vizMap = zeros(rows, cols);
    vizMap(isHotPixel) = 1;      % Mark all hot
    vizMap(isIsolated) = 0.5;    % Mark isolated distinct from clustered
    
    imagesc(vizMap);
    % Custom colormap: Black (Normal), Green (Isolated), Red (Clustered)
    colormap(gca, [0 0 0; 0 1 0; 1 0 0]); 
    % title(sprintf('Hot Pixel Clusters (Red=Clustered, Green=Isolated)'), 'FontSize', mainFontSize);
    
    title(sprintf('Red=Clustered, Green=Isolated'), 'FontSize', mainFontSize);
    xlabel('Sensor X', 'FontSize', mainFontSize);
    set(gca, 'FontSize', mainFontSize);
    axis image;

    exportgraphics(f3, fullfile('img', 'dark_hot_pixels.png'), 'Resolution', 300);

    % --- FIGURE 3: Hot Pixel Map (VISIBILITY FIX) ---
    f4 = figure('Name', 'Hot Pixel Map', 'Color', 'w', 'Units', 'pixels', ...
        'Position', [200 200 figWidth figHeight]);
    
    hold on;
    % Plot black background 
    imagesc(ones(rows, cols));
    colormap(gca, 'gray'); 
    axis image; set(gca, 'YDir', 'reverse'); % Ensure correct image orientation
    
    % 2. Find coordinates for scatter plot
    [y_iso, x_iso] = find(isIsolated);
    [y_clst, x_clst] = find(isHotPixel & ~isIsolated); % Clustered are hot but not isolated
    
    if ~isempty(x_iso)
        scatter(x_iso, y_iso, 30, [0 1 0], 'filled'); % Green for Isolated
    end
    if ~isempty(x_clst)
        scatter(x_clst, y_clst, 30, [1 0 0], 'filled'); % Red for Clustered
    end
    
    title(sprintf('Red=Clustered, Green=Isolated'), 'FontSize', mainFontSize);
    xlabel('Sensor X', 'FontSize', mainFontSize);
    ylabel('Sensor Y', 'FontSize', mainFontSize);
    set(gca, 'FontSize', mainFontSize);
    
    % Force limits to match sensor size
    xlim([1 cols]);
    ylim([1 rows]);
    hold off;

   exportgraphics(f4, fullfile('img', 'dark_hot_pixels_big.png'), 'Resolution', 300);

    
end
end