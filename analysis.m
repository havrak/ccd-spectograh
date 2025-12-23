function [imageStack] = getImageStack(folder)
files = dir(fullfile(folder, '*.fit'));
N = numel(files);
imageStack = [];


for i = 1:N
	image = files(i);
	data = fitsread(fullfile(image.folder, image.name)); % we don't have to crop the images
	imageStack(:,:,i) = double(data);
end
end


function [imageStack] = getImageStackExposure(folder, exposure)
files = dir(fullfile(folder, '*.fit'));
N = numel(files);
imageStack = [];

index = 1;
for i = 1:N
	image = files(i);
	metadata = fitsinfo(fullfile(image.folder, image.name));


	[metaExposure, ~] = scripts.getMetadataRow(metadata.PrimaryData.Keywords, 'EXPTIME');

	if(metaExposure ~= exposure)
		continue
	end
	data = fitsread(fullfile(image.folder, image.name)); % we don't have to crop the images
	imageStack(:,:,index) = double(data);
	index = index+1;
end
end

clc;

zeroStack = getImageStack("data/zero/"); % they have two different exposures (one is one second longer) but the difference is minimal
zeroResults = analyzeZeros(zeroStack);
displayZero = 0;
displayDarks = 0;

if displayZero
figure(1);
subplot(2,2,1);
imagesc(zeroResults.masterZero);
colorbar; title('Master Bias Frame');

subplot(2,2,2);
imagesc(zeroResults.readNoiseMap);
colorbar; title('Read Noise Map');

subplot(2,2,3);
histogram(zeroResults.masterZero(:), 100);
title('Bias Level Distribution');
xlabel('ADU'); ylabel('Count');

subplot(2,2,4);
histogram(zeroResults.readNoiseMap(:), 100);
title('Read Noise Distribution');
xlabel('ADU'); ylabel('Count');

fprintf('Global Read Noise: %.6f ADU\n', zeroResults.readNoise);
end




darkStack = getImageStack("data/dark/");
exposureTimeDarks = 3600; % exact time doesn't matter
darkResults = analyzeDarks(darkStack, zeroResults, exposureTimeDarks);


if displayDarks
    figure(2);
    subplot(2,2,1);
    imagesc(darkResults.masterDark);
    colorbar; title('Master Dark Frame (median)');
    
    %subplot(2,2,2);
    %imagesc(darkResults.darkCurrentMap * 3600); % electrons/hour
    %colorbar; title('Dark Current (e-/hour)');
    
    subplot(2,2,2);
    imagesc(darkResults.darkNoiseMap);
    colorbar; title('Read + Dark Current Noise');
    
    subplot(2,2,3);
    histogram(darkResults.darkCurrentMap(:) * 3600, 100);
    title('Dark Current Distribution');
    xlabel('e-/hour'); ylabel('Count');
    
    subplot(2,2,4);
    histogram(darkResults.darkNoiseMap(:), 50);
    title('Read + Dark Current Noise');
    xlabel('ADU'); ylabel('Count');
    
    % subplot(2,2,4);
    % histogram(darkResults.shotNoiseMap(:), 50);
    % title('Shot Noise Distribution (removed readNoise var)');
    % xlabel('ADU'); ylabel('Count');
    
    fprintf('Mean Dark Current: %.6f e-/s/pixel\n', darkResults.meanDarkCurrent);
    fprintf('Mean Dark Noise: %.2f ADU\n', darkResults.meanDarkNoise);
end


% analyzeFlatfieldsInterFrame("data/flatfield_whole/", zeroResults, darkResults);
% analyzeFlatfieldsSpatial("data/flatfield_one_conf/", zeroResults, darkResults);
% analyzeFlatfieldsSpatialRobust("data/flatfield_one_conf/", zeroResults, darkResults);
% analyzeFlatfieldsInterFrameRobust("data/flatfield_whole/", zeroResults, darkResults);

% analyzeFlatfieldsModel("data/flatfield_one_conf/", zeroResults);
% analyzeFlatfieldsInterFrameModel("data/flatfield_one_conf/", zeroResults);

analyzeFlatfieldsPixelModel("data/flatfield_whole/", zeroResults);