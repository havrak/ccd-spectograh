function zeroResults = analyzeZeros(folder, plotFlag)
zeroStack = scripts.getImageStack(folder);

if nargin < 2
    plotFlag = 0;
end

masterZero = median(zeroStack, 3); % create a reference zero frame

readNoiseMap = var(double(zeroStack), 0, 3);
readNoiseMapStd = sqrt(readNoiseMap); 

validNoiseMap = readNoiseMap(readNoiseMap > 1e-6);

globalReadNoiseVar = median(validNoiseMap); 
globalReadNoiseVal = sqrt(globalReadNoiseVar);

mu_bias = mean(masterZero(:));
sigma_bias = std(double(masterZero(:)));

zeroResults.masterZero = masterZero;              % master zero frame
zeroResults.meanLevel = mu_bias;                  % mean of zero frame  
zeroResults.meanLevelStd = sigma_bias;            % what is deviation within zero frame
zeroResults.readNoiseMap = readNoiseMapStd;       % what is deviation within all zeros
zeroResults.globalReadNoise = globalReadNoiseVal; % meadian standard deviation of read noise
 
fprintf('------------------------------------------------\n');
fprintf('ZERO FRAMES ANALYSIS RESULTS\n');
fprintf('------------------------------------------------\n');
fprintf('Processed:             %d frames\n', size(zeroStack, 3));
fprintf('Mean Bias Level:       %.4f ADU\n', mu_bias);
fprintf('Bias Level Std Dev:    %.4f ADU (Spatial Non-uniformity)\n', sigma_bias); 
fprintf('Global Read Noise:     %.4f ADU (Temporal RMS)\n', globalReadNoiseVal);
fprintf('------------------------------------------------\n\n');

if plotFlag
    figure('Name', 'Zero Frame Analysis', 'Color', 'w');

    subplot(2,2,1);
    imagesc(zeroResults.masterZero);
    colorbar;
    title('Master Bias Frame (Median)');
    %axis image; 
    colormap('gray');

    subplot(2,2,2);
    imagesc(zeroResults.readNoiseMap);
    colorbar;
    title('Read Noise Map (RMS)');
    %axis image;
    clim([0, prctile(zeroResults.readNoiseMap(:), 99)]);


    subplot(2,2,3);
    hold on;
    data_bias = double(masterZero(:));
    histogram(data_bias, 'Normalization', 'pdf', 'EdgeColor', 'none');

    x_fit = linspace(min(data_bias), max(data_bias), 100);
    y_fit = normpdf(x_fit, mu_bias, sigma_bias);

    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    title(sprintf('Bias Uniformity (\\mu=%.2f, \\sigma=%.2f)', mu_bias, sigma_bias));
    xlabel('ADU'); ylabel('Probability Density');
    legend('Data', 'Gaussian Fit');
    hold off;

    subplot(2,2,4);
    valid_noise_vector = readNoiseMapStd(readNoiseMapStd > 0.1);

    histogram(valid_noise_vector, 50, 'EdgeColor', 'none');
    title('Read Noise Distribution (Zeros Excluded)');
    xlabel('Read Noise (ADU RMS)'); ylabel('Count');

    xline(globalReadNoiseVal, 'r--', 'LineWidth', 1.5, 'Label', 'Median');
end
end