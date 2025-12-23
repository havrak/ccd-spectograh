function calibrateSpectraBulk(dataFolder, thFile, arFile)
% CalibrateSpectraBulk - Global calibration using angle constraints.
% 1. Loads all headers and spectra.
% 2. Separates by Spectral Order.
% 3. Optimizes Dispersion and Angle-Wavelength relationship globally.

% ---------------------------------------------------------
% STAGE 1: LOAD REFERENCE SPECTRA
% ---------------------------------------------------------
fprintf('Loading reference data...\n');
[refLam, refInt] = loadAndMerge(thFile, arFile);

% Define the Master Grid (0.05 Angstrom resolution)
% This is the "Common Resolution" we convert everything to for comparison.
dLambda = 0.05;
gridAxis = 3000:dLambda:10000;
refSpectrum = zeros(size(gridAxis));

% Gaussian Smoothing
for i = 1:length(refLam)
    if refLam(i) < 3000 || refLam(i) > 10000, continue; end
    idx = round((refLam(i) - 3000) / dLambda) + 1;
    refSpectrum(idx) = refSpectrum(idx) + refInt(i);
end
refSpectrum = smoothdata(refSpectrum, 'gaussian', round(0.5/dLambda));
refSpectrum = refSpectrum / max(refSpectrum);

% ---------------------------------------------------------
% STAGE 2: LOAD IMAGES
% ---------------------------------------------------------
files = dir(fullfile(dataFolder, '*.fit*'));
dataset = [];

fprintf('Loading %d files...\n', length(files));
for k = 1:length(files)
    fname = files(k).name;
    [~, meta, data] = scripts.getMetadata(dataFolder, fname);
    fprintf("Image %s, Gtating %4.2f, Order %2.0f, Exposure %3.0f, Type %s\n", fname, meta.gratingAng, meta.order, meta.exposure, meta.type);

    if meta.type ~= "comp"
        continue;
    end

    % Create 1D Spectrum (95th percentile reduces cosmic ray impact)
    spec1D = prctile(data, 95, 1);
    spec1D = spec1D - min(spec1D); % Baseline removal
    spec1D = spec1D / max(spec1D); % Normalize

    entry.name = fname;
    entry.spectrum = spec1D;
    entry.angle = meta.gratingAng;
    entry.pos = meta.gratingPos;
    entry.order = meta.order;
    entry.exp = meta.exposure;
    entry.nx = length(spec1D);

    % Basic Peak Detection for the Voting Algorithm
    [pks, locs] = findpeaks(spec1D, 'MinPeakProminence', 0.05, 'MinPeakDistance', 5);
    [~, pIdx] = sort(pks, 'descend');
    nKeep = min(40, length(pks));
    entry.peaksVals = pks(pIdx(1:nKeep));
    entry.peaksPix = locs(pIdx(1:nKeep));

    dataset = [dataset; entry];

end

% ---------------------------------------------------------
% STAGE 3: PROCESS IMAGES IN SPECTRAL ORDER
% ---------------------------------------------------------
runOrderSearch(dataset, 1, refSpectrum, gridAxis, dLambda, refLam, refInt);
runOrderSearch(dataset, 2, refSpectrum, gridAxis, dLambda, refLam, refInt);
end

function runOrderSearch(dataset, order, refSpec, gridAxis, dLam, refRawLam, refRawInt)
idx = [dataset.order] == order;
subset = dataset(idx);

if isempty(subset)
    fprintf('\nNo images found for Spectral Order %d.\n', order);
    return;
end

fprintf('Processing Order %d (%d images)\n', order, length(subset));


% ---------------------------------------------------------
% STAGE 2: DEFINE A SEARCH PARAMETERS
% ---------------------------------------------------------
if order == 1 % Red
    slopeEst = (8635 - 5565) / (35.78 - 27.88); % ~388 A/deg
    dispRange = 0.15:0.005:0.80;
    slopeRange = slopeEst * (0.95:0.01:1.05);
else % Blue
    slopeEst = (4874 - 4390) / (38.75 - 36.17); % ~187 A/deg
    dispRange = 0.10:0.001:0.50; % [0.1142578125];
    slopeRange = slopeEst * (0.95:0.01:1.05);
end

% ---------------------------------------------------------
% STAGE 3: PERFORM THE SEARCH
% ---------------------------------------------------------
results = [];
gridMin = min(gridAxis);
fprintf('Grid Search: %d Slopes x %d Disps...\n', length(slopeRange), length(dispRange));

% Pre-computation for speed
nGrid = length(gridAxis);

for s = slopeRange
    for d = dispRange

        % --- STACKING STEP ---
        % Construct "Hypothesis Spectrum" for this Slope/Dispersion
        hypSpec = zeros(1, nGrid);

        for k = 1:length(subset)
            img = subset(k);

            % Calculate Relative Wavelengths (Assuming Intercept 0)
            % CenterWav = Angle * Slope
            % PixelWav  = CenterWav + Disp * (Pix - Center)
            pix = 1:img.nx;
            relWav = (s * img.angle) + d * (pix - img.nx/2);

            % Map to Grid Indices
            gridIdx = round((relWav - gridMin) / dLam) + 1;
            valid = gridIdx >= 1 & gridIdx <= nGrid;

            % Accumulate
            hypSpec(gridIdx(valid)) = hypSpec(gridIdx(valid)) + img.spectrum(valid);
        end

        % --- CORRELATION STEP ---
        if max(hypSpec) == 0, continue; end
        hypSpec = hypSpec / max(hypSpec);

        % xcorr finds alignment (lags)
        [c, lags] = xcorr(refSpec, hypSpec);

        [pks, peakIdx] = findpeaks(c, 'SortStr', 'descend', 'NPeaks', 5, 'MinPeakDistance', 20);

        % Store all found peaks for this geometry
        for p = 1:length(pks)
            thisLag = lags(peakIdx(p));
            interceptFound = thisLag * dLam;
            results = [results; pks(p), s, d, interceptFound];
        end
    end
end

% ---------------------------------------------------------
% STAGE 4: SELECT TOP 5
% ---------------------------------------------------------
[~, sortIdx] = sort(results(:,1), 'descend');
top = results(sortIdx(1:min(10, end)), :);

fprintf('\n--- TOP 5 CANDIDATES ---\n');
fprintf('  # | Score  | Slope (A/deg) | Disp (A/px) |  Reso (A) | Intercept (A)\n');
fprintf('----------------------------------------------------------------------\n');
for i = 1:size(top, 1)
    fprintf('  %d | %.4f | %13.4f | %11.4f | %9.4f | %13.1f\n', ...
        i, top(i,1), top(i,2), top(i,3), top(i,3)*2048, top(i,4));
end

% best = top(2,:);
% visualizeResult(subset, best(2), best(3), best(4), refRawLam, refRawInt);
end

function visualizeResult(subset, slope, disp, intercept, refLam, refInt)
figure('Name', sprintf('Order Result: Slope=%.2f Disp=%.3f', slope, disp));

subplot(2,1,2);
stem(refLam, refInt/max(refInt), 'r', 'Marker', 'none');
xlabel('Wavelength (A)'); ylabel('Reference Intensity');
title('ThAr Reference'); grid on;
hold on;

subplot(2,1,1);
hold on;
colors = jet(length(subset));

minW = inf; maxW = -inf;

fprintf('\n--- INDIVIDUAL IMAGE STATS (Best Match) ---\n');
fprintf('  %-25s | Angle | Center (A)\n', 'Image Name');
fprintf('  -------------------------------------------------\n');

for k = 1:length(subset)
    img = subset(k);
    centerLam = (slope * img.angle) + intercept;
    pix = 1:img.nx;
    lamAxis = centerLam + disp * (pix - img.nx/2);

    plot(lamAxis, img.spectrum, 'Color', colors(k,:), 'LineWidth', 1);
    minW = min(minW, min(lamAxis));
    maxW = max(maxW, max(lamAxis));


    fprintf('  %-25s | %5.2f | %10.2f f\n', ...
        img.name, img.angle, centerLam);
end

title(sprintf('Stacked Observations (Slope: %.2f, Disp: %.3f)', slope, disp));
ylabel('Normalized Intensity'); grid on;
xlim([minW, maxW]);

linkaxes([subplot(2,1,1), subplot(2,1,2)], 'x');
subplot(2,1,2); xlim([minW, maxW]);
end

function [mergedLam, mergedInt] = loadAndMerge(thPath, arPath)
[thInt, thLam] = cleanParse(thPath);
[arInt, arLam] = cleanParse(arPath);
allInt = [thInt; arInt];
allLam = [thLam; arLam];
[mergedLam, idx] = sort(allLam);
mergedInt = allInt(idx);
[mergedLam, uIdx] = unique(mergedLam);
mergedInt = mergedInt(uIdx);
end

function [intensity, lambda] = cleanParse(filePath)
str = fileread(filePath);
str = regexprep(str, '\', '');
str = regexprep(str, '[a-zA-Z]+', ' ');
nums = sscanf(str, '%f');
if mod(length(nums), 2) ~= 0, nums = nums(1:end-1); end
data = reshape(nums, 2, []).';
intensity = data(:, 1);
lambda = data(:, 2);
end



% first for 5100A - 9100A and
% second for 3700A - 5100A
calibrateSpectraBulk("data/thar_all","data/thoriumtable2_a.txt", "data/argontable2_a.txt");
