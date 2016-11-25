function [Fmaps] = fillMapofTwoFaceShapes( name1, name2, numEigBasis )
numEig = numEigBasis;

% lazy loading, till error computation
mesh1 = loadMeshLB( name1, numEig);
mesh2 = loadMeshLB( name2, numEig);

segMethod.name = 'manu';
segMethod.faceSeg = 0;

%% use mouth boundray
s=load([name1, '_segs.mat']);
%ml = 2*ones(mesh1.nv, 1);
%ml(mouth, :) = 1;
mesh1.segmentation = s.seg;
mesh1.segMethod = segMethod;

s=load([name2, '_segs.mat']);
%ml = 2*ones(mesh2.nv, 1);
%ml(mouth, :) = 1;
mesh2.segmentation = s.seg;
mesh2.segMethod = segMethod;

%%%%%% Basis
% use conf LB basis - unweighted cot matrix
basis1 = mesh1.confLaplaceBasis(:,1:numEig);
basis2 = mesh2.confLaplaceBasis(:,1:numEig);

%%%%%% Descriptors
% descriptor options
options.descriptors = {};
options.descriptorsParams = {};

% Point descriptors
% WKS
options.descriptors{end+1}          = 'wks';
wksParams.timeSteps = 100;
wksParams.type = 'shape';
wksParams.number = wksParams.timeSteps;
options.descriptorsParams{end+1}    = wksParams;      

% % WKM
options.descriptors{end+1}          = 'wkm';
wkmParams.timeSteps = 10;
wkmParams.type = 'shape';
wkmParams.number = 10;
match1 = load([name1, '_', num2str(wkmParams.number), 'pt.mat']);
wkmParams.landmarks(:,1) = match1.match;
match2 = load([name2, '_', num2str(wkmParams.number), 'pt.mat']);
wkmParams.landmarks(:,2) = match2.match;
options.descriptorsParams{end+1}    = wkmParams;   


segMatches = [[1],[1]];    
options.descriptors{end+1}          = 'segment_indicator';
segParams.segs(:,1) = [1 2 3 4 5 6 7];
segParams.segs(:,2) = [1 2 3 4 5 6 7];
options.descriptorsParams{end+1}    = segParams; 

% Texture
% TODO
 options.descriptors{end+1}          = 'rgb';
 colourParams.type = 'texture';
 colourParams.number = 3;
 options.descriptorsParams{end+1}    = colourParams; 

% Compute descriptors
fprintf(1, 'Computing descriptors...\n');
[f, g] = collectDescriptors(mesh1, mesh2, options);

% Write descriptors in provided basis
fprintf(1, 'Writing descriptors in provided basis...\n');
F = descriptorsToCoeffs(f, basis1);
G = descriptorsToCoeffs(g, basis2);

fprintf(1, 'Re-weighting descriptors...\n');
[F, G] = reweighDescriptors(F, G, numEig, options);

%%%%%%%%% Operator commutativity
options.operators = {};

% LB operator
options.operators{end+1} = 'confLB';
%options.operators{end+1} = 'LB';

% collect operators
fprintf(1, 'Computing operators...\n');
[s, r] = collectOperators(mesh1, mesh2, options);

% Write operators in provided basis
fprintf(1, 'Writing operators in provided basis...\n');
S = operatorsToMtx(s, basis1);
R = operatorsToMtx(r, basis2);

for m=1:length(R)
    if (nnz(S{m}) > 10*mesh1.nv || nnz(R{m}) > 10*mesh2.nv)
        error('?');
    end
end

%%%%%%%%% Optimization
% Which maps. options are: band, groundTruth, l1, leastSquares, qr, svd
options.mappingMethods = {};
options.mappingParams = {};

options.mappingMethods{end+1} = 'leastSquares';
options.mappingParams{end+1} = {};

% Perform ICP in eigenvec space for 
options.postprocessing = 'ICP';
options.postprocessingParams = {10};               % step:10


Fmaps{1,1} = [];
Fmaps{1,1}.method = options.mappingMethods;
Fmaps{1,1}.params = options.mappingParams;
Fmaps{1,1}.ppMethod = options.postprocessing;
Fmaps{1,1}.ppParams = options.postprocessingParams;
Fmaps{1,1}.numEig =  numEig;
Fmaps{1,1}.segMatches = segMatches;
Fmaps{1,1}.segMethod = segMethod;
Fmaps{1,1}.F = F;
Fmaps{1,1}.G = G;
Fmaps{1,1}.S = S{1};
Fmaps{1,1}.R = R{1};
Fmaps{1,1}.mesh1 = [];
Fmaps{1,1}.mesh2 = [];
Fmaps{1,1}.basis1 = basis1;
Fmaps{1,1}.basis2 = basis2;
Fmaps{1,1}.options = options;

clear mesh1 mesh2
