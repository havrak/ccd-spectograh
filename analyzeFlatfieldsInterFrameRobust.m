function gainResults = analyzeFlatfieldsInterFrameRobust(folder, zeroResults, darkResults)

% outlier rejection parameters
numBins = 50;
outlierSigma = 2.0;

files = dir(fullfile(folder, '*.fit'));
rnSquared = zeroResults.globalReadNoise^2;

% ---------------------------------------------------------
% STAGE 1: GROUPING
% ---------------------------------------------------------
groupSettings = containers.Map('KeyType', 'char', 'ValueType', 'any');

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
fprintf('Calculating statistics from image pairs...\n');
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

	mask =  mu < 55000;
	validMeans = mu(mask);
	validVars  = sigma2(mask);

	allMeans = [allMeans; validMeans];
	allVars  = [allVars; validVars];
end

% ---------------------------------------------------------
% STAGE 3: BINNING & ROBUST REJECTION
% ---------------------------------------------------------

fprintf('Binning %d raw points and rejecting outliers...\n', length(allMeans));

binnedSignal = zeros(numBins, 1);
binnedVar    = zeros(numBins, 1);
binnedStd    = zeros(numBins, 1);

edges = linspace(0, max(allMeans), numBins+1);

for i = 1:numBins
	idxs = allMeans >= edges(i) & allMeans < edges(i+1);

	if sum(idxs) < 10, continue; end

	binVars = allVars(idxs);
	binMeans = allMeans(idxs);

	medVar = median(binVars);
	madVar = median(abs(binVars - medVar));
	sigmaEst = 1.4826 * madVar; % Convert MAD to Sigma

	upper = medVar + (outlierSigma * sigmaEst);
	lower = medVar - (outlierSigma * sigmaEst);

	cleanMask = (binVars <= upper) & (binVars >= lower);

	binnedSignal(i) = mean(binMeans(cleanMask));
	binnedVar(i)    = mean(binVars(cleanMask));
	binnedStd(i)    = std(binVars(cleanMask));
end

% Remove empty bins
validBins = binnedSignal > 0;
xData = binnedSignal(validBins);
yData = binnedVar(validBins);
yStd  = binnedStd(validBins);

% ---------------------------------------------------------
% STAGE 4: CALCULATE GAIN
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
% STAGE 5: VISUALIZATION (THE SQUARE ROOT PLOT)
% ---------------------------------------------------------

figure(101); clf;
allNoise = sqrt(allVars);
scatter(allMeans(1:50:end), allNoise(1:50:end), 2, 'k', 'filled', 'MarkerFaceAlpha', 0.05);
hold on;

errorbar(xData, sqrt(yData), sqrt(yStd)./sqrt(sqrt(length(allMeans)/numBins)), 'ko');

xFit = linspace(0, max(xData), 100);
yFit = sqrt( (1/gainVal)*xFit + rnSquared );
plot(xFit, yFit, 'r-', 'LineWidth', 2);

title(sprintf('Photon Transfer Curve (Inter-Frame Pair Difference)\nGain = %.3f e-/ADU', gainVal));
xlabel('Mean Signal (ADU)');
ylabel('Noise (RMS ADU)');
grid on;
axis tight;;
xlim([0 60000]);
legend('Binned Data', sprintf('Fit (G=%.2f)', gainVal), 'Location', 'NorthWest');

% ---------------------------------------------------------
% RESULTS
% ---------------------------------------------------------
gainResults.gain = gainVal;
gainResults.readNoiseSquared = rnSquared;

fprintf('------------------------------------------------\n');
fprintf('INTER-FRAME ANALYSIS RESULTS (ROBUST)\n');
fprintf('------------------------------------------------\n');
fprintf('Gain Estimated:        %.4f e-/ADU\n', gainVal);
fprintf('Fixed Intercept (RN^2): %.2f ADU^2\n', rnSquared);
fprintf('Calc Intercept (RN^2): %.2f ADU^2\n', intercept);
fprintf('------------------------------------------------\n');
end
