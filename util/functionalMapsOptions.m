function options = functionalMapsOptions()
% Creates a default functional mapping options struct.
%%%%%%%%%%%%% not used!

options = [];
options.timeSteps = 30;                     % hks/wks timesteps
options.descriptors = {'wks','hkm','wkm'};  % which descriptors

% % landmarks
% options.landmarks1 = landmarks;
% options.landmarks2 = landmarks;

% The user can specify files with additional descriptors in
% options.descriptorFiles (one for each mesh)


% Which maps. options are: band, groundTruth, l1, l2, leastSquares, qr, svd
options.mappingMethods = {};
options.mappingParams = {};
% L1, param = lambda. see l1RegularizedMap
options.mappingMethods{end+1} = 'l1';
options.mappingParams{end+1}  = {1};
% SVD, param = epsilon, see svdMap
options.mappingMethods{end+1} = 'svd';
options.mappingParams{end+1}  = {1e-4};
