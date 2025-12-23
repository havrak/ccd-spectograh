function gainResults = analyzeFlatfieldsSpatial(folder, zeroResults, darkResults)

% binning paramters
binY = 5;

% stratified sampling paramters
numBins = 50;
samplesPerBin = 2000;

files = dir(fullfile(folder, '*.fit'));
rnSquared = zeroResults.globalReadNoise^2;

% ---------------------------------------------------------
% STAGE 1: SPATIAL VARIANCE CALCULATION
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
% STAGE 2: STRATIFIED SAMPLING
% ---------------------------------------------------------
% Instead of random sampling, we bin the data by intensity and
% take equal samples from each bin. This forces high-signal regions
% to be included in the fit.

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
% STAGE 3: CALCULATE GAIN
% ---------------------------------------------------------
% Math: Variance = (1/Gain) * Mean + ReadNoise^2
fprintf('Fitting Gain Curve to %d spatial points...\n', length(allMeans));

Y_fit = allMeans - rnSquared;
X_fit = allVars;

b = robustfit(X_fit, Y_fit, 'bisquare', 4.685);
slope = b(2);
intercept = b(1);

gainVal = 1 / slope;

% ---------------------------------------------------------
% STAGE 4: VISUALIZATION
% ---------------------------------------------------------
figure(101); clf;

allNoise = sqrt(allVars);
scatter(allMeans, allNoise, 2, 'k', 'filled', 'MarkerFaceAlpha', 0.05);
hold on;

xFit = linspace(0, 60000, 100);
yFit = sqrt( (1/gainVal) * xFit + intercept );
plot(xFit, yFit, 'r-', 'LineWidth', 2);

title(sprintf('Photon Transfer Curve (Spatial Method)\nGain = %.3f e-/ADU, bin size: X=1 Y=%.3f', gainVal, binY));
xlabel('Mean Signal (ADU)');
ylabel('Noise (RMS ADU)');
grid on;
axis tight;
xlim([0 60000]);
legend('Measured Data', ['Fit: \sigma = \surd(S / ' sprintf('%.2f', gainVal) ')']);


% ---------------------------------------------------------
% RESULTS
% ---------------------------------------------------------
gainResults.gain = gainVal;
gainResults.readNoiseSquared = intercept;

fprintf('------------------------------------------------\n');
fprintf('SPATIAL ANALYSIS RESULTS\n');
fprintf('------------------------------------------------\n');
fprintf('Gain Estimated:        %.4f e-/ADU\n', gainVal);
fprintf('Read Noise (fit):      %.2f ADU\n', sqrt(abs(intercept)));
fprintf('------------------------------------------------\n');
end
