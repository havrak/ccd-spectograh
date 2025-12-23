function darkAnalysis = analyzeDarks(darkStack, zeroResults, exposureTime)
	% Subtract bias from darks
	%   variance of pixels will still be effected by both readNoise and
	%   shotNoise, but we remove a basic pattern from it
	%
	darksCorrected = darkStack - zeroResults.masterZero;


	dataMedian = median(darksCorrected(:));
	threshold = 2000;

	darksCorrected(darksCorrected > threshold) = dataMedian;


	% Calculate master dark
	masterDark = median(darksCorrected, 3);

	% Dark current calculation (per second)
	darkCurrent = masterDark / exposureTime; % electron/ADU generation per second

	% Noise analysis
	darkNoise = std(darksCorrected, 0, 3);

	% Separate shot noise from read noise
	% Total variance = read_noise + shot_noise

	%readNoiseVariance = zeroResults.globalReadNoise^2;
	%totalVariance = var(darksCorrected, 0, 3);

	darkAnalysis.masterDark = masterDark;
	darkAnalysis.darkCurrentMap = darkCurrent;
	darkAnalysis.darkNoiseMap = darkNoise;
	%darkAnalysis.shotNoiseMap = sqrt(max(shotNoiseVariance, 0)); % avoid negative values, corrected wit

	% Global statistics
	darkAnalysis.meanDarkCurrent = mean(darkCurrent(:));
	darkAnalysis.meanDarkNoise = mean(darkNoise(:));

	% Plot results

end
