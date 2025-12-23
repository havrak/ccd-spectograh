function [data, metadata, dataUncropped] = getMetadata(folder, image)

dataUncropped = fitsread(fullfile(folder, image));
metadataWhole = fitsinfo(fullfile(folder, image));
metadata = [];

[metadata.exposure, ~] = scripts.getMetadataRow(metadataWhole.PrimaryData.Keywords, 'EXPTIME');


[specFilt, ~] = scripts.getMetadataRow(metadataWhole.PrimaryData.Keywords, 'SPECFILT');

if specFilt == 5
    metadata.order = 2; % 3700A - 5100A 
elseif specFilt == 1
    metadata.order = 1; % 5100A - 9100A 
else
    metadata.order = 0xFF;
end

[metadata.gratingPos, ~] = scripts.getMetadataRow(metadataWhole.PrimaryData.Keywords, 'GRATPOS');

[metadata.gratingAng, ~] = scripts.getMetadataRow(metadataWhole.PrimaryData.Keywords, 'GRATANG');

[metadata.type, ~] = scripts.getMetadataRow(metadataWhole.PrimaryData.Keywords, 'OBJECT');

[crop, ~] = scripts.getMetadataRow(metadataWhole.PrimaryData.Keywords, 'TRIMSEC');

crop = crop(2:end-1);

parts = strsplit(crop, ',');

xRange = strsplit(parts{1}, ':');
x1 = str2double(xRange{1});
x2 = str2double(xRange{2});

yRange = strsplit(parts{2}, ':');
y1 = str2double(yRange{1});
y2 = str2double(yRange{2});

metadata.crop.x1 = x1;
metadata.crop.x2 = x2;

metadata.crop.y1 = y1;
metadata.crop.y2 = y2;

data = dataUncropped(y1:y2, x1:x2);

% figure(1);
% subplot(212);
% imagesc(dataUncropped);
% subplot(211);
% imagesc(data);
end


