function zeroResults = analyzeZeros(zeroStack)
	% Calculate master zero/bias frame
	masterZero = median(zeroStack, 3);

	% Analyze read noise
	readNoise = var(zeroStack, 0, 3);

	% Statistical analysis
	zeroResults.masterZero = masterZero;
	zeroResults.readNoiseMap = readNoise;
	zeroResults.readNoise = mean(readNoise(:));
	zeroResults.readNoiseStd = sqrt(readNoise(:)); % std from std???
end

