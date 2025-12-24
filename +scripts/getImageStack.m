function [imageStack] = getImageStack(folder)
files = dir(fullfile(folder, '*.fit'));
N = numel(files);
imageStack = [];


for i = 1:N
	image = files(i);
	data = fitsread(fullfile(image.folder, image.name)); 
	imageStack(:,:,i) = double(data);
end
end