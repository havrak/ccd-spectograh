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
hotPixelThreshold = meanDark + (10 * stdDarkSpatial);
isHotPixel = darkCurrentMap_e_s > hotPixelThreshold;
numHotPixels = sum(isHotPixel(:));

% Identify Transient Events (Gamma/Beta Rays)
% Gamma Ray is a transient event (High in ONE frame -> High Std, Low Median)
maxFrame = max(darksCorrected, [], 3);

isCosmicRay = (maxFrame - masterDark) > transientThreshold;
numCosmicRays = sum(isCosmicRay(:));

% Transient events would mess with statistics, we will just sacrifice them
% and replace their values with median
darksCleaned = darksCorrected;
for i = 1:numFrames
    thisFrame = darksCorrected(:,:,i);

    diffMap = abs(thisFrame - masterDark);
    isTransient = diffMap > transientThreshold;

    thisFrame(isTransient) = masterDark(isTransient);
    darksCleaned(:,:,i) = thisFrame;
end


% Noise Analysis
observedNoise = std(darksCorrected, 0, 3);
observedNoiseCleaned = std(darksCleaned, 0, 3);

% Results
darkAnalysis.masterDark = masterDark;
darkAnalysis.darkCurrentMap = darkCurrentMap_e_s;   % in e-
darkAnalysis.observedNoiseMap = observedNoise;    % in e-
darkAnalysis.observedNoiseCleanedMap = observedNoiseCleaned;

darkAnalysis.hotPixelCount = numHotPixels;
darkAnalysis.meanDarkCurrent = mean(darkCurrentMap_e_s(:));
darkAnalysis.meanDarkNoise = mean(observedNoise(:));
darkAnalysis.meanDarkNoiseCleaned = mean(observedNoiseCleaned(:));

fprintf('------------------------------------------------\n');
fprintf('DARK FRAMES ANALYSIS RESULTS (16-bit)\n');
fprintf('------------------------------------------------\n');
fprintf('Analysis Configuration:\n');
fprintf('  Frames:            %d\n', numFrames);
fprintf('  Exposure Time:     %.1f s\n', exposureTime);
fprintf('  Gain:              %.2f e-/ADU\n', gain);
fprintf('------------------------------------------------\n');
fprintf('Master Dark Current Statistics:\n');
fprintf('  Mean Dark Noise:   %.7f e- (RMS)\n', stdDarkSpatial);
fprintf('  Mean Dark Noise:   %.2f ADU (RMS)\n', std(masterDark(:)));
fprintf('------------------------------------------------\n');
fprintf('Sensor Statistics:\n');
fprintf('  Mean Dark Current: %.4f e-/s/pixel\n', darkAnalysis.meanDarkCurrent);
fprintf('  Mean Dark Noise:   %.2f ADU (RMS)\n', darkAnalysis.meanDarkNoise);
fprintf('  Mean Dark Noise:   %.2f ADU (RMS - Cleaned)\n', darkAnalysis.meanDarkNoiseCleaned);
fprintf('------------------------------------------------\n');
fprintf('Defects & Radiation:\n');
fprintf('  Static Hot Pixels: %d (%.3f%% of sensor)\n', numHotPixels, (numHotPixels/(rows*cols))*100);
fprintf('  Gamma/Cosmic Hits: %d (Transient Events)\n', numCosmicRays);
fprintf('------------------------------------------------\n\n');

if plotFlag
    figure('Name', 'Dark Frame Analysis', 'Color', 'w');

    % 1. Master Dark (Cleaned)
    subplot(1,3,1);
    imagesc(darkAnalysis.masterDark);
    colorbar;
    colormap('gray');
    title('Master Dark (Gamma Removed)');
    clim([0 prctile(masterDark(:), 99)]);

    % 2. Dark Current Hist (Log Scale)
    subplot(1,3,2);
    histogram(darkAnalysis.darkCurrentMap(:), 100);
    set(gca, 'YScale', 'log');
    title('Dark Current Distribution');
    xlabel('Dark Current (e-/s)'); ylabel('Count (Log)');
    grid on;

    % 3. Hot Pixel Map
    subplot(1,3,3);
    imagesc(isHotPixel);
    colormap(gca, 'gray');
    title(sprintf('Static Hot Pixels (> %.2f e-/s)', hotPixelThreshold));
    xlabel('White = Defect');
end
end