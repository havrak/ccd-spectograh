function gainResults = analyzeFlatfieldsPixelModel(folder, zeroResults, darkResults)
% Estimates Gain using a "Residuals" approach:
% 1. Build a "Model" frame for each fiber configuration using the average of the most common exposure.
% 2. For every frame, scale the Model to match the frame's global mean.
% 3. Calculate residuals (Actual - Model).
% 4. Analyze Variance of residuals vs Signal to determine Gain.

BIN_SETTINGS.numBins      = 1000;   % Bins for the final PTC plot
BIN_SETTINGS.minPoints    = 100;    % Min points to trust a final bin
BIN_SETTINGS.outlierSigma = 2.5;   % Outlier rejection strictness
BIN_SETTINGS.maxADU       =45000; % Saturation cutoff

% Sampling Settings (Per Image)
SAMPLE.numLevels     = 500;  % How many intensity levels to split each image into
SAMPLE.pointsPerLevel = 10000; % Max points to take per level (20k points max per image
rnSquared = zeroResults.globalReadNoise.^2; % Use pre-calculated Read Noise

INTEROLATION_SETTING.calcReadNoise = 0;
INTEROLATION_SETTING.statisticsToolbox = 1;

files = dir(fullfile(folder, '*.fit'));
if isempty(files)
    error('No FITS files found in folder: %s', folder);
end

% ---------------------------------------------------------
% STAGE 1: PARSE AND GROUP FILES
% ---------------------------------------------------------
fprintf('Parsing file structure...\n');
fileData = struct();
fibreMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:numel(files)
    fileData(i).name = files(i).name;
    fileData(i).path = fullfile(folder, files(i).name);

    parts = split(files(i).name, '_');
    if length(parts) >= 3
        fileData(i).exposure = str2double(parts{2});
        fileData(i).fibre = parts{3};

        if isKey(fibreMap, fileData(i).fibre)
            fibreMap(fileData(i).fibre) = [fibreMap(fileData(i).fibre), i];
        else
            fibreMap(fileData(i).fibre) = [i];
        end
    end
end

fibKeys = keys(fibreMap);
fprintf('Found %d fiber configurations.\n', length(fibKeys));

% ---------------------------------------------------------
% STAGE 2: BUILD MODELS & CALCULATE RESIDUALS
% ---------------------------------------------------------
fprintf('Building pixel models and calculating residuals...\n');

allSignals = [];    % Expected Signal (from Model)
allResiduals = [];  % (Actual - Model)
globalMeans = [];   % For Plot 2
globalExps = [];    % For Plot 2
globalFibreIdx = [];% For Plot 2 coloring

% Helper to load image
loadImg = @(path) double(fitsread(path)) - zeroResults.masterZero;

for f = 1:length(fibKeys)
    key = fibKeys{f};
    indices = fibreMap(key);

    % --- Step 2a: Find "Reference" Exposure for this Fibre ---
    % We want the exposure time that has the MOST frames to average out noise.
    % expTimes = [fileData(indices).exposure];
    % uniqueExps = unique(expTimes);
    % counts = histcounts(expTimes, [uniqueExps, max(uniqueExps)+1]);
    %
    % [~, bestIdx] = max(counts);
    % refExpTime = uniqueExps(bestIdx);
    %
    % % Find indices of these reference frames
    % refIndices = indices(expTimes == refExpTime);
    %
    % fprintf('  Fibre %s: Using %d frames at %.2fs as reference model.\n', ...
    %     key, length(refIndices), refExpTime);
    %
    % % -- Step 2b: Create Master Reference Frame ---
    % refStack = [];
    % for k = 1:length(refIndices)
    %     refStack(:,:,k) = loadImg(fileData(refIndices(k)).path);
    % end


    % --- Step 2a: Create Master Reference Frame ---
    refStack = [];
    for k = 1:length(indices) % use all frames to calculate the mean
        refStack(:,:,k) = loadImg(fileData(indices(k)).path);
    end


    % Calculate the Pattern (normalized model)
    masterModel = double(mean(refStack, 3));
    masterModel(masterModel < zeroResults.globalReadNoise*5 ) = 0; % nuke non illuminated pixels that would than get scaled

    meanRef = mean(masterModel(:));
    normalizedPattern = masterModel ./ meanRef;

    % --- Step 2b: Calculate Residuals for ALL frames in this group ---
    for k = 1:length(indices)
        idx = indices(k);
        img = double(loadImg(fileData(idx).path));

        % NEW: Regress the Image against the Pattern (Fit the Scale Factor)
        % We only fit pixels that are part of the signal (masking the dark background)
        % This prevents background noise/bias from skewing the scaling factor.

        fitMask = normalizedPattern > zeroResults.globalReadNoise*5/ meanRef; % Use only illuminated features

        if sum(fitMask(:)) > 100
            % Solve linear equation: img = scaleFactor * pattern
            % The '\' operator performs robust least-squares regression
         ratios = img(fitMask) ./ normalizedPattern(fitMask);
             
             % The Scaling Factor is the Median of these ratios.
             % This is immune to background drifts and outliers.
             scaleFactor = median(ratios);
        else
            % Fallback for empty/dark frames
            scaleFactor = mean(img(:));
        end

        % Create Prediction using the fitted scalar
        prediction = normalizedPattern .* scaleFactor;

        % 1. Determine Scale Factor based on MEAN INTENSITY
        currentMean = mean(img(:));

        % Store for Plot 2
        globalMeans = [globalMeans; currentMean];
        globalExps  = [globalExps; fileData(idx).exposure];
        globalFibreIdx = [globalFibreIdx; f];


        % 3. Calculate Residual
        residual = img - prediction;

        % 4. Select valid pixels
        % - Ignore saturation
        % - Ignore very dark pixels (dead pixels)
        % - Ignore edges (optional, usually handled by crop in getMetadata but we do full frame here)
        mask = prediction > 100 & prediction < BIN_SETTINGS.maxADU;

        validPred = prediction(mask);
        validRes  = residual(mask);

        if isempty(validPred)
            continue;
        end

        allSignals = [allSignals; validPred];
        allResiduals = [allResiduals; validRes];
    end
end

[allSignals, allResiduals] = stratifiedSample(allSignals, allResiduals, SAMPLE.numLevels, SAMPLE.pointsPerLevel);


% ---------------------------------------------------------
% STAGE 3: BINNING & STATISTICS (PTC GENERATION)
% ---------------------------------------------------------
fprintf('Analyzing %d points...\n', numel(allSignals));

[allSignals, sortIdx] = sort(allSignals);
allResiduals = allResiduals(sortIdx);

binEdges = linspace(min(allSignals), max(allSignals), BIN_SETTINGS.numBins+1);
binnedSignals = [];
binnedVars  = []; % This will hold Mean Squared Residuals
binnedStds  = []; % RMS of Residuals

for b = 1:BIN_SETTINGS.numBins
    % 1. Isolate data in this signal range
    inBin = allSignals >= binEdges(b) & allSignals < binEdges(b+1);
    sigs = allSignals(inBin);
    resids = allResiduals(inBin);

    if length(sigs) < BIN_SETTINGS.minPoints
        continue;
    end

    % 2. Robust Outlier Rejection (MAD)
    % We filter based on deviation from the bin median
    medRes = median(resids);
    madRes = median(abs(resids - medRes));
    sigmaEst = 1.4826 * madRes;

    validMask = abs(resids - medRes) < (BIN_SETTINGS.outlierSigma * sigmaEst);
    cleanResids = resids(validMask);
    cleanSigs   = sigs(validMask);

    %if isempty(cleanResids), continue; end

    % 3. Calculate "Real" Deviance Statistics
    % We do NOT use var() or std() because they subtract the local mean.
    % We assume the Model is true, so residual IS the noise.

    meanSqResid = mean(cleanResids.^2); % Pure Variance (E[x^2])
    rmsResid    = sqrt(meanSqResid);    % Pure Noise (RMS)

    binnedSignals = [binnedSignals; mean(cleanSigs)];
    binnedVars  = [binnedVars; meanSqResid];
    binnedStds  = [binnedStds; rmsResid];
end

fprintf("Fitting into %.0f points\n", length(binnedSignals));

% ---------------------------------------------------------
% STAGE 4: FITTING GAIN
% ---------------------------------------------------------

% --- Fit per bin ---
% Equation: gain*Std = sqrt((1/Gain) * Signal + (1/Gain^2) * ReadNoise^2
% y = mx + c
% m = 1/Gain
% c = ReadNoise^2 (We can fix this or let it float. Let's fix it for robustness)

% var = signal/gain + readNoise
localGain = binnedSignals ./ (binnedVars - rnSquared);

% --- Polyfit of whole function ---
fprintf('Fitting full noise model...\n');
X = binnedSignals;
Y = binnedVars;

g0 = 1.0;
rn0 = zeroResults.globalReadNoise;
options = optimset('Display','off', 'TolX', 1e-6);

if INTEROLATION_SETTING.statisticsToolbox == 1
    if INTEROLATION_SETTING.calcReadNoise == 1
        modelFun = @(p, x) x ./ p(1) + (p(2)^2)./(p(1).^2);
        [beta, R, J, CovB, MSE] = nlinfit(X, Y, modelFun, [g0, rn0]);
        ci = nlparci(beta, R, 'Jacobian', J);

        % Standard Errors (Half-width of CI / 1.96 approx, or sqrt(diag(CovB)))
        se = sqrt(diag(CovB));
        gainVal = beta(1);
        gainSE  = se(1);  % Standard Error of Gain
        fitRN   = beta(2);
        rnSE    = se(2);  % Standard Error of RN
        p_opt = [gainVal, fitRN];
    else
        modelFun = @(p, x) x ./ p(1) + rn0^2/p(1)^2;
        [beta, R, J, CovB, MSE] = nlinfit(X, Y, modelFun, [g0]);
        ci = nlparci(beta, R, 'Jacobian', J);

        % Standard Errors (Half-width of CI / 1.96 approx, or sqrt(diag(CovB)))
        se = sqrt(diag(CovB));
        gainVal = beta(1);
        gainSE  = se(1);  % Standard Error of Gain
        fitRN   = zeroResults.globalReadNoise;
        p_opt = [gainVal];
        rnSE    = 0;  % Standard Error of RN
    end

    fprintf('------------------------------------------------\n');
    fprintf('FLATFIELD ANALYSIS RESULTS\n');
    fprintf('------------------------------------------------\n');
    fprintf('Gain value:            %.6f +/- %.4f e-/ADU\n', gainVal, 1.96*gainSE);
    fprintf('Noise value:           %.4f +/- %.2f ADU\n', fitRN, 1.96*rnSE);
    fprintf('------------------------------------------------\n\n');

else
    if INTEROLATION_SETTING.calcReadNoise == 1
        modelFun = @(p, x) x ./ p(1) + (p(2)^2)./(p(1).^2);
        objective = @(p) sum( (Y - modelFun(p, X)).^2 );
        p_opt = fminsearch(objective, [g0, rn0], options);

        fitRN   = p_opt(2);
    else
        modelFun = @(p, x) x ./ p(1) + rn0/p(1)^2;
        objective = @(p) sum( (Y - modelFun(p, X)).^2 );
        p_opt = fminsearch(objective, [g0], options);
        fitRN   = zeroResults.globalReadNoise;

    end

    gainVal = p_opt(1);
    fprintf('------------------------------------------------\n');
    fprintf('FLATFIELD ANALYSIS RESULTS\n');
    fprintf('------------------------------------------------\n');
    fprintf('Gain value:            %.6f e-/ADU\n', gainVal);
    fprintf('Noise value:           %.4f ADU\n', fitRN);
    fprintf('------------------------------------------------\n\n');
end



% ---------------------------------------------------------
% STAGE 5: SHUTTER OFFSET ESTIMATION
% ---------------------------------------------------------
% We estimate this from the Mean vs Exposure plot data
% Model: Signal = Flux * (Time + Offset)
% Time = Signal/Flux - Offset
% We do a simple robust fit over all fibers
shutterOffsets = [];
for f = 1:length(fibKeys)
    fIdx = globalFibreIdx == f;
    fMeans = globalMeans(fIdx);
    fExps  = globalExps(fIdx);

    if length(unique(fExps)) > 1
        p = polyfit(fExps, fMeans, 1); %
        % x-intercept is when y=0 -> 0 = p1*x + p2 -> x = -p2/p1
        % This x is -(shutter_offset).
        % So Offset = p2/p1

        % p = robustfit(fExps, fMeans);
        offset = p(2) / p(1);
        shutterOffsets = [shutterOffsets; offset];
    end
end
avgShutterOffset = median(shutterOffsets);

% ---------------------------------------------------------
% STAGE 6: VISUALIZATION
% ---------------------------------------------------------

figWidth = 1000;
figHeight = 800;
mainFontSize = 18;



f1 = figure('Name', 'Master Bias', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [100 100 figWidth figHeight]);
scatter(binnedSignals, binnedStds, 15, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
xFit = linspace(0, max(binnedSignals), 200);
yFit = sqrt(modelFun(p_opt, xFit));
plot(xFit, yFit, 'r-', 'LineWidth', 2);
plot(xFit, sqrt(xFit ./ gainVal), 'b:', 'LineWidth', 1.5);
xlabel('Signal (ADU)'); ylabel('Noise \sigma (ADU)');
%title('Photon Transfer Curve');
set(gca, 'FontSize', mainFontSize);
grid on;
legend('Data', sprintf('Fit (G=%.4f)', gainVal), 'Shot Noise Only', 'Location', 'southeast');
exportgraphics(f1, fullfile('img', 'model_photon_transfer_curve.png'), 'Resolution', 300);

figure(); clf;
t = tiledlayout(2, 2, 'TileSpacing', 'compact');
title(t, sprintf('Gain Analysis (Residual Method) | Est. Gain: %.3f e-/ADU', gainVal));

% --- PLOT 0: Noise Fit ---
nexttile;
scatter(binnedSignals, binnedStds, 15, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
xFit = linspace(0, max(binnedSignals), 200);
yFit = sqrt(modelFun(p_opt, xFit));
plot(xFit, yFit, 'r-', 'LineWidth', 2);
plot(xFit, sqrt(xFit ./ gainVal), 'b:', 'LineWidth', 1.5);
xlabel('Signal (ADU)'); ylabel('Noise \sigma (ADU)');
title('Photon Transfer Curve');
grid on; legend('Data', sprintf('Fit (G=%.4f)', gainVal), 'Shot Noise Only', 'Location', 'southeast');


% --- PLOT 1: Photon Transfer Curve ---
nexttile;
scatter(binnedSignals, binnedVars, 10, 'k', 'filled');
hold on;
% Plot Fit
xFit = linspace(0, max(binnedSignals), 100);
yFit = (1/gainVal) * xFit + rnSquared;
plot(xFit, yFit, 'r-', 'LineWidth', 1.5);
xlabel('Signal (ADU)');
ylabel('Variance (ADU^2)');
title('Photon Transfer Curve');
grid on;
legend('Binned Data', 'Gain Fit');

% --- PLOT 2: Gain vs Signal (Linearity) ---
nexttile;
plot(binnedSignals, localGain, 'b.-');
hold on;
yline(gainVal, 'r--', 'LineWidth', 2);
xlabel('Signal Intensity (ADU)');
ylabel('Instantaneous Gain (e-/ADU)');
title('Gain Linearity');
grid on;
ylim([0, gainVal*2]); % Zoom in relevant area
legend('Calculated Gain per Bin', 'Global Fit Gain');

% --- PLOT 3: Mean Value vs Exposure Time ---
nexttile;
hold on;
colors = lines(length(fibKeys));

for f = 1:length(fibKeys)
    % Extract data for this fiber
    fIdx = globalFibreIdx == f;
    fExps = globalExps(fIdx);
    fMeans = globalMeans(fIdx);

    % 1. Plot RAW data as discrete dots (no connecting lines)
    scatter(fExps, fMeans, 20, colors(f,:), 'filled', 'MarkerFaceAlpha', 0.6);

    % 2. Calculate average mean for each unique exposure time
    uniExps = unique(fExps);
    avgMeans = zeros(size(uniExps));
    for u = 1:length(uniExps)
        avgMeans(u) = mean(fMeans(fExps == uniExps(u)));
    end

    % 3. Plot the TREND line through the averages
    [sExp, sI] = sort(uniExps);
    plot(sExp, avgMeans(sI), '-', 'Color', colors(f,:), 'LineWidth', 1.5);
end

xlabel('Exposure Time (s)');
ylabel('Mean Signal (ADU)');
title(sprintf('Linearity Check (Offset \\approx %.3fs)', avgShutterOffset));
grid on;

% --- PLOT 4: Residual Distribution (Histogram) ---

% fprintf("All residuals: ");
% disp(sliceResids);

f2 = figure('Name', 'Master Bias', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [100 100 figWidth figHeight]);

medSignal = median(allSignals);

disp(medSignal);
sliceWindow = 2000; % ADU width
sliceMask = abs(allSignals - medSignal) < sliceWindow;
sliceResids = allResiduals(sliceMask);


h = histogram(sliceResids, 300, 'Normalization', 'pdf');
hold on;

% Overlay Theoretical Gaussian
% Width should be sqrt(Signal/Gain + RN^2)
theoreticalSigma = sqrt( (medSignal / gainVal) + rnSquared );
xH = linspace(min(sliceResids), max(sliceResids), 100);
disp(size(xH));
yH = normpdf(xH, 0, theoreticalSigma);

plot(xH, yH, 'r-', 'LineWidth', 2);
xlabel('Residual (ADU)');
ylabel('Probability Density');
title(sprintf('@ ~%d ADU (Exp \\sigma=%.1f)', round(medSignal), theoreticalSigma));
legend('Measured Residuals', 'Poisson Model');
grid on;
xlim([-4*theoreticalSigma, 4*theoreticalSigma]);

set(gca, 'FontSize', mainFontSize);
% exportgraphics(f2, fullfile('img', 'model_histogram_med.png'), 'Resolution', 300);


f3 = figure('Name', 'Master Bias', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [100 100 figWidth figHeight]);

q = quantile(allSignals, [0.25 0.5 0.75]); % 1st, 2nd (median), 3rd quartiles
medSignal = q(3);

disp(medSignal);
sliceWindow = 2000; % ADU width
sliceMask = abs(allSignals - medSignal) < sliceWindow;
sliceResids = allResiduals(sliceMask);

% fprintf("All residuals: ");
% disp(sliceResids);

h = histogram(sliceResids, 300, 'Normalization', 'pdf');
hold on;

% Overlay Theoretical Gaussian
% Width should be sqrt(Signal/Gain + RN^2)
theoreticalSigma = sqrt( (medSignal / gainVal) + rnSquared );
xH = linspace(min(sliceResids), max(sliceResids), 100);
disp(size(xH));
yH = normpdf(xH, 0, theoreticalSigma);

plot(xH, yH, 'r-', 'LineWidth', 2);
xlabel('Residual (ADU)');
ylabel('Probability Density');
title(sprintf('@ ~%d ADU (Exp \\sigma=%.1f)', round(medSignal), theoreticalSigma));
legend('Measured Residuals', 'Poisson Model');
grid on;
xlim([-4*theoreticalSigma, 4*theoreticalSigma]);


set(gca, 'FontSize', mainFontSize);
% exportgraphics(f3, fullfile('img', 'model_histogram_3qar.png'), 'Resolution', 300);


f4 = figure('Name', 'Master Bias', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [100 100 figWidth figHeight]);
colors = lines(length(fibKeys));
hold on;

% 1. Create an array to store the specific handles for the legend
scatPlots = gobjects(1, length(fibKeys)); 

for f = 1:length(fibKeys)
    % Extract data for this fiber
    fIdx = globalFibreIdx == f;
    fExps = globalExps(fIdx);
    fMeans = globalMeans(fIdx);
    
    % 2. Assign the output of scatter to the array
    scatPlots(f) = scatter(fExps, fMeans, 20, colors(f,:), 'filled', 'MarkerFaceAlpha', 0.6);
    
    % Calculate average mean for each unique exposure time
    uniExps = unique(fExps);
    avgMeans = zeros(size(uniExps));
    for u = 1:length(uniExps)
        avgMeans(u) = mean(fMeans(fExps == uniExps(u)));
    end
    
    % 3. Plot the TREND line (we don't save this handle, so legend ignores it)
    [sExp, sI] = sort(uniExps);
    plot(sExp, avgMeans(sI), '-', 'Color', colors(f,:), 'LineWidth', 1.5);
end

xlabel('Exposure Time (s)');
ylabel('Mean Signal (ADU)');
title(sprintf('Offset \\approx %.3fs', avgShutterOffset));

% 4. Pass the specific 'scatPlots' handles to the legend function
legend(scatPlots, ['Intensity A'; 'Intensity B'; 'Intensity C'; 'Intensity D'], ...
    'Location','northeast');

grid on;

set(gca, 'FontSize', mainFontSize);
exportgraphics(f4, fullfile('img', 'exposure_time.png'), 'Resolution', 300);

f5 = figure('Name', 'Master Bias', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [100 100 figWidth figHeight]);
 plot(binnedSignals, localGain, 'b.-');
hold on;
yline(gainVal, 'r--', 'LineWidth', 2);
xlabel('Signal Intensity (ADU)');
ylabel('Instantaneous Gain (e-/ADU)');
grid on;
ylim([0, gainVal*2]); % Zoom in relevant area
legend('Calculated Gain per Bin', 'Global Fit Gain');

set(gca, 'FontSize', mainFontSize);
% exportgraphics(f5, fullfile('img', 'model_instantaneous_gain.png'), 'Resolution', 300);

% ---------------------------------------------------------
% RETURN
% ---------------------------------------------------------
gainResults.gain = gainVal;
gainResults.shutterOffset = avgShutterOffset;
gainResults.readNoiseSquared = rnSquared;

end


function [sampledSig, sampledRes] = stratifiedSample(signals, residuals, numBins, maxPerBin)

if isempty(signals)
    sampledSig = [];
    sampledRes = [];
    return;
end

sampledSig = [];
sampledRes = [];

minVal = min(signals);
maxVal = max(signals);

binEdges = linspace(minVal, maxVal, numBins + 1);

[~, ~, binIndices] = histcounts(signals, binEdges);

for b = 1:numBins
    idxs = find(binIndices == b);

    if isempty(idxs)
        continue;
    end

    count = length(idxs);
    if count > maxPerBin
        perm = randperm(count, maxPerBin);
        keepIdx = idxs(perm);
    else
        keepIdx = idxs;
    end

    sampledSig = [sampledSig; signals(keepIdx)];
    sampledRes = [sampledRes; residuals(keepIdx)];
end
end