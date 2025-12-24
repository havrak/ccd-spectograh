function gainResults = analyzeFlatfieldsInterFrame(folder, zeroResults)
files = dir(fullfile(folder, '*.fit'));

% stratified sampling paramters
numBins = 50;
samplesPerBin = 2000;

% ---------------------------------------------------------
% STAGE 1: GROUPING
% ---------------------------------------------------------
groupSettings = containers.Map('KeyType', 'char', 'ValueType', 'any');
groupFibre    = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:numel(files)
	nameSplit = split(files(i).name, '_');
	expStr = nameSplit{2};
	fibStr = nameSplit{3};

	keySet = sprintf('%s_%s', expStr, fibStr);
	if isKey(groupSettings, keySet)
		groupSettings(keySet) = [groupSettings(keySet), i];
	else
		groupSettings(keySet) = [i];
	end
end

% ---------------------------------------------------------
% STAGE 2: CALCULATE RAW STATISTICS (Pairwise Difference)
% ---------------------------------------------------------
fprintf('Calculating pixel statistics...\n');
allMeans = [];
allVars  = [];

keys = groupSettings.keys;
for k = 1:length(keys)
	indices = groupSettings(keys{k});

	% We need at least 2 images to verify noise
	if length(indices) < 2
		fprintf("Discarding entry: %s\n", keys{k});
		continue;
	end

	% Load Stack
	stack = [];
	for j = 1:length(indices)
		idx = indices(j);
		data = double(fitsread(fullfile(files(idx).folder, files(idx).name)));
		data = data - zeroResults.masterZero;
		stack(:,:,j) = data;
	end

	% Calculate Pixel Statistics
	mu = mean(stack, 3);
	sigma2 = var(stack, 0, 3); % instead of std we use var so that final equation will be linear

	mask = mu > 0 & mu < 60000;
	validMeans = mu(mask);
	validVars  = sigma2(mask);

	allMeans = [allMeans; validMeans];
	allVars  = [allVars; validVars];
end

% ---------------------------------------------------------
% STAGE 3: STRATIFIED SAMPLING
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
% STAGE 4: CALCULATE GAIN
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
% STAGE 5: VISUALIZATION (THE SQUARE ROOT PLOT)
% ---------------------------------------------------------
figure(101); clf;
allNoise = sqrt(allVars);
scatter(allMeans, allNoise, 1, 'k', 'filled', 'MarkerFaceAlpha', 0.1);
hold on;

xFit = linspace(min(allMeans), max(allMeans), 100);
yFit = sqrt( (1/gainVal) * xFit + intercept );
plot(xFit, yFit, 'r-', 'LineWidth', 2);

title(sprintf('Photon Transfer Curve (Noise vs Signal)\nGain = %.3f e-/ADU', gainVal));
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
gainResults.shutterOffset = avgOffset;
gainResults.readNoiseSquared = intercept;

fprintf('------------------------------------------------\n');
fprintf('INTER-FRAME ANALYSIS RESULTS\n');
fprintf('------------------------------------------------\n');
fprintf('Gain Estimated:        %.4f e-/ADU\n', gainVal);
fprintf('Fit Intercept (RN^2):  %.2f\n', intercept);
fprintf('------------------------------------------------\n');
end
