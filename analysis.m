clc;

zeroResults = analyzeZeros("data/zero/",0);
darkResults = analyzeDarks("data/dark/", zeroResults, 0.623635, 1);



% analyzeFlatfieldsInterFrame("data/flatfield_whole/", zeroResults, darkResults);
% analyzeFlatfieldsSpatial("data/flatfield_one_conf/", zeroResults, darkResults);
% analyzeFlatfieldsSpatialRobust("data/flatfield_one_conf/", zeroResults, darkResults);
% analyzeFlatfieldsInterFrameRobust("data/flatfield_whole/", zeroResults, darkResults);

% analyzeFlatfieldsModel("data/flatfield_one_conf/", zeroResults);
% analyzeFlatfieldsInterFrameModel("data/flatfield_one_conf/", zeroResults);

% analyzeFlatfieldsPixelModel("data/flatfield_whole/", zeroResults);