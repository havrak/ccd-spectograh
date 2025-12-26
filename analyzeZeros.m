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
    % Check/Create Export Directory
    if ~exist('img', 'dir')
        mkdir('img');
    end

    % Configuration for consistent export look
    % We use a square-ish aspect ratio which fits well in document columns
    figWidth = 1000;
    figHeight = 800;
    mainFontSize = 18; % Large font to remain readable when image is resized small

    % --- FIGURE 1: Master Bias ---
    f1 = figure('Name', 'Master Bias', 'Color', 'w', 'Units', 'pixels', ...
        'Position', [100 100 figWidth figHeight]);

    imagesc(zeroResults.masterZero);
    colorbar;
    % title('Master Bias Frame (Median)', 'FontSize', mainFontSize);
    colormap(gca, 'gray');
    set(gca, 'FontSize', mainFontSize);
    axis image; % Keep aspect ratio of the sensor
   xlabel('Sensor X', 'FontSize', mainFontSize);
    ylabel('Sensor Y', 'FontSize', mainFontSize);
    exportgraphics(f1, fullfile('img', 'zero_master_frame.png'), 'Resolution', 300);

    % --- FIGURE 2: Read Noise Map ---
    f2 = figure('Name', 'Read Noise Map', 'Color', 'w', 'Units', 'pixels', ...
        'Position', [150 150 figWidth figHeight]);

    imagesc(zeroResults.readNoiseMap);
    colorbar;
    % title('Read Noise Map (RMS)', 'FontSize', mainFontSize);
    colormap(gca, 'gray');

    clim([0, prctile(zeroResults.readNoiseMap(:), 99)]); % Smart contrast
    set(gca, 'FontSize', mainFontSize);
    axis image;
       xlabel('Sensor X', 'FontSize', mainFontSize);
    ylabel('Sensor Y', 'FontSize', mainFontSize);

    exportgraphics(f2, fullfile('img', 'zero_noise_map.png'), 'Resolution', 300);

    % --- FIGURE 3: Master Bias Histogram ---
    f3 = figure('Name', 'Master Bias Hist', 'Color', 'w', 'Units', 'pixels', ...
        'Position', [200 200 figWidth figHeight]);

    hold on;
    data_bias = double(masterZero(:));
    histogram(data_bias,'Normalization', 'pdf', 'EdgeColor', 'none');

    x_fit = linspace(min(data_bias), max(data_bias), 100);
    y_fit = normpdf(x_fit, mu_bias, sigma_bias);

    plot(x_fit, y_fit, 'r-', 'LineWidth', 3); % Thicker line for visibility
    title(sprintf('\\mu=%.2f, \\sigma=%.2f', mu_bias, sigma_bias), 'FontSize', mainFontSize);
    % title(sprintf('Master Uniformity\n\\mu=%.2f, \\sigma=%.2f', mu_bias, sigma_bias), 'FontSize', mainFontSize);
    xlabel('ADU', 'FontSize', mainFontSize);
    ylabel('Probability Density', 'FontSize', mainFontSize);
    legend({'Data', 'Gaussian Fit'}, 'FontSize', mainFontSize, 'Location', 'best');
    set(gca, 'FontSize', mainFontSize);
    grid on;
    hold off;

    exportgraphics(f3, fullfile('img', 'zero_master_uniformity.png'), 'Resolution', 300);

    % --- FIGURE 4: Full Stack Histogram ---
    f4 = figure('Name', 'Stack Noise Hist', 'Color', 'w', 'Units', 'pixels', ...
        'Position', [250 250 figWidth figHeight]);

    hold on;
    allPixels = double(zeroStack(:));

    % Outlier rejection for statistics (avoids skewing the fit)
    mask_outliers = allPixels > 625;
    allPixels(mask_outliers) = median(allPixels(~mask_outliers));

    mu_all = mean(allPixels);
    sigma_all = std(allPixels);

    histogram(allPixels, 'BinWidth', 1, 'Normalization', 'pdf', ...
        'EdgeColor', 'none', 'FaceColor', [0.8500 0.3250 0.0980]);

    x_fit_stack = linspace(min(allPixels), max(allPixels), 100);
    y_fit_stack = normpdf(x_fit_stack, mu_all, sigma_all);

    plot(x_fit_stack, y_fit_stack, 'k--', 'LineWidth', 3);

    title(sprintf('\\mu=%.2f, \\sigma=%.2f', mu_all, sigma_all), 'FontSize', mainFontSize);
    % title(sprintf('Total Stack Noise\n\\mu=%.2f, \\sigma=%.2f', mu_all, sigma_all), 'FontSize', mainFontSize);
    xlabel('ADU', 'FontSize', mainFontSize);
    ylabel('Probability Density', 'FontSize', mainFontSize);
    legend({'All Pixels', 'Gaussian Fit'}, 'FontSize', mainFontSize, 'Location', 'best');
    set(gca, 'FontSize', mainFontSize);
    grid on;
    hold off;

    exportgraphics(f4, fullfile('img', 'zero_stack_uniformity.png'), 'Resolution', 300);
end
end