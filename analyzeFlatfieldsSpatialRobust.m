function gainResults = analyzeFlatfieldsSpatialRobust(folder, zeroResults, darkResults)

% binning paramters
binY = 5;

% outlier rejection
outlierSigma  = 5.0;
numSignalBins = 50;

% stratified sampling paramters
numBins = 50;
samplesPerBin = 1000;

files = dir(fullfile(folder, '*.fit'));
rnSquared = zeroResults.globalReadNoise^2;

% ---------------------------------------------------------
% STAGE 1: RAW SPATIAL VARIANCE
% ---------------------------------------------------------
fprintf('Calculating Spatial Variance (Single Image Method)...\n');
fprintf('Binning: X=1 Y=%.3f\n', binY);
fprintf('Analyzing %d files\n', numel(files));

allMeans = [];
allVars  = [];

dofFactor = (binY - 1) / (binY - 2); % correction for reduced degrees of freedom

for i = 1:numel(files)
	data = double(fitsread(fullfile(folder, files(i).name)));
	data = data-zeroResults.masterZero;

	[rows, cols] = size(data);

	% Vectorized Sliding Window Analysis
	% We process column-by-column for speed
	for col = 1:cols
		colData = data(:, col);

		% Truncate to multiple of binY
		nChunks = floor(rows / binY);
		if nChunks < 1, continue; end

		% Reshape: Each column of 'chunks' is a spatial group
		chunks = reshape(colData(1:nChunks*binY), binY, nChunks);

		% Detrend: Remove vertical intensity gradients (light profile)
		% chunksDetrended = detrend(chunks);
		% varVal = var(chunksDetrended).*dofFactor; % Noise Level

		% Statistics
		mu = mean(chunks); % detrend shoudn't affect mu
		varVal = var(chunks).*dofFactor;

		mask = mu > 0 & mu < 60000;
		allMeans = [allMeans; mu(valid)'];
		allVars  = [allVars; varVal(valid)'];
	end
end

% ---------------------------------------------------------
% STAGE 2: BINNING & CLEANING
% ---------------------------------------------------------
fprintf('Performing binning and cleaning...\n');


binnedSignal = zeros(numSignalBins, 1);
binnedVar    = zeros(numSignalBins, 1);
binnedStd    = zeros(numSignalBins, 1);

edges = linspace(0, max(allMeans), numSignalBins+1);

for k = 1:numSignalBins
	idxs = allMeans >= edges(k) & allMeans < edges(k+1);
	if sum(idxs) < 5, continue; end

	binVars = allVars(idxs);
	binMeans = allMeans(idxs);

	medVar = median(binVars);
	madVar = median(abs(binVars - medVar));
	sigmaEst = 1.4826 * madVar;

	cutoffUpper = medVar + outlierSigma * sigmaEst;
	cutoffLower = medVar - outlierSigma * sigmaEst;

	cleanIdx = (binVars <= cutoffUpper) & (binVars >= cutoffLower);

	binnedSignal(k) = mean(binMeans(cleanIdx));
	binnedVar(k)    = mean(binVars(cleanIdx));
	binnedStd(k)    = std(binVars(cleanIdx));
end

% Remove empty bins
validBins = binnedSignal > 0;
xData = binnedSignal(validBins);
yData = binnedVar(validBins);
yStd  = binnedStd(validBins);

% ---------------------------------------------------------
% STAGE 3: CALCULATE GAIN
% ---------------------------------------------------------
% Math: Variance = (1/Gain) * Mean + ReadNoise^2
fprintf('Fitting Gain Curve to %d spatial points...\n', length(xData));

Y_fit = yData - rnSquared;
X_fit = xData;

b = robustfit(X_fit, Y_fit, 'bisquare', 4.685);
slope = b(2);
intercept = b(1);

gainVal = 1 / slope;

% ---------------------------------------------------------
% STAGE 4: STRATIFIED SAMPLING
% ---------------------------------------------------------
% Here this is just to get a managable amount of data for the visualization

fprintf('Performing Stratified Sampling...\n');

finalMeans = [];
finalVars = [];

binEdges = linspace(min(allMeans), max(allMeans), numBins);
[~, ~, binIndices] = histcounts(allMeans, binEdges);

for b = 1:length(binEdges)-1
	idxs = find(binIndices == b);

	if isempty(idxs)
		continue;
	end

	% If we have more points than needed, sample randomly within this bin
	if length(idxs) > samplesPerBin
		perm = randperm(length(idxs), samplesPerBin);
		selectedIdxs = idxs(perm);
	else
		selectedIdxs = idxs;
	end

	finalMeans = [finalMeans; allMeans(selectedIdxs)];
	finalVars  = [finalVars; allVars(selectedIdxs)];
end

allMeans = finalMeans;
allVars = finalVars;

% ---------------------------------------------------------
% STAGE 5: VISUALIZATION
% ---------------------------------------------------------
figure(101); clf;

allNoise = sqrt(allVars);
scatter(allMeans, allNoise, 2, 'k', 'filled', 'MarkerFaceAlpha', 0.05);
hold on;

errorbar(xData, sqrt(yData), sqrt(yStd)./sqrt(sqrt(length(allMeans)/numSignalBins)), 'ko', 'MarkerFaceColor', 'k');

xFit = linspace(0, max(xData), 200);
yFit = sqrt( (1/gainVal)*xFit + rnSquared );
plot(xFit, yFit, 'r-', 'LineWidth', 2);

title(sprintf('Photon Transfer Curve (Noise Scale)\nGain = %.3f e-/ADU', gainVal));
xlabel('Mean Signal (ADU)');
ylabel('Noise (RMS ADU)');
grid on;
axis tight;
xlim([0 60000]);
legend('Measured Data','Binned Data', sprintf('Fit (G=%.2f)', gainVal), 'Fit Range Limit', 'Location', 'NorthWest');


figure(102); clf;
if ~isempty(sampleBinVars)
	histogram(sampleBinVars, 50, 'Normalization', 'pdf'); hold on;
	xline(sampleBinCutoff, 'r-', 'LineWidth', 2, 'Label', sprintf('Sigma Limit (%.1f)', outlierSigma));

	title('DIAGNOSTIC: Variance Distribution (Sample Bin)');
	xlabel('Pixel Variance'); ylabel('Probability Density');
end

% ---------------------------------------------------------
% RESULTS
% ---------------------------------------------------------
gainResults.gain = gainVal;
gainResults.x = xData;
gainResults.y_noise = yDataNoise;

fprintf('------------------------------------------------\n');
fprintf('SPATIAL ANALYSIS RESULTS\n');
fprintf('------------------------------------------------\n');
fprintf('Gain Estimated:        %.4f e-/ADU\n', gainVal);
fprintf('Read Noise (fit):      %.2f ADU\n', sqrt(abs(intercept)));
fprintf('------------------------------------------------\n');

end
