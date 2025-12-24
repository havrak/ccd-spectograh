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
