function Fmaps  = computeMap(basis1, F, S, basis2, G, R, options)
% Computes a functional map between using basis1 and basis2, the
% descriptors F and G, and the operators in S and R 
% Note that some mapping methods ignore the operators 

Fmaps = [];
for i=1:length(options.mappingMethods)
    method = options.mappingMethods{i};
    params = options.mappingParams{i};
    
    fprintf(1, 'Map: %s ...', method);
    t = tic;

    Fmaps{i} = doMap(basis1, F, S, basis2, G, R, method, params);
    
    fprintf(1, '\t took: %g seconds\n', toc(t));
end

% Do postprocessing
if isfield(options,'postprocessing')
    for i=1:length(Fmaps)
        pp = options.postprocessing;
        ppParams = options.postprocessingParams;
        Fmaps{i} = doPostprocessing(Fmaps{i}, basis1, basis2, pp, ppParams);
    end
end

% Compute maps using each of the specified methods
function Fmap = doMap(basis1, F, S, basis2, G, R, method, params)
Fmap.method = method;
Fmap.params = params;

% Band. Params: band
% NOTE: doesn't use operators
if strcmp(method,'band')
    for j=1:length(params)
        Fmap.maps{j} = bandDiagonalMap(F,G,params{j});
    end

% Ground truth. Params: none
elseif strcmp(method,'groundTruth'),
    Fmap.maps{1} = groundTruthMap(basis1,basis2);


% L1. Params: lambda
% NOTE: doesn't use operators
elseif strcmp(method,'l1')
    for j=1:length(params)
        Fmap.maps{j} = l1RegularizedMap(F,G,params{j});
    end

% L2. Params: lambda
% NOTE: doesn't use operators
elseif strcmp(method,'l2')
    for j=1:length(params)
        Fmap.maps{j} = l2RegularizedMap(F,G,params{j});
    end

% Least squares. Params: none
elseif strcmp(method,'leastSquares'),
    Fmap.maps{1} = leastSquaresMapFast(F,G,S,R);

% QR. Params: epsilon    
% NOTE: doesn't use operators
elseif strcmp(method,'qr'), 
    for j=1:length(params)
        Fmap.maps{j} = qrMap(F,G, params{j});
    end

% SVD. Params: epsilon    
% NOTE: doesn't use operators
elseif strcmp(method,'svd')
    for j=1:length(params)
        Fmap.maps{j} = svdMap(F,G,params{j});
    end

else
    error(['Unrecognized mapping method: ' method]);
end

% Do some postprocessing
function outmaps = doPostprocessing(inmaps, basis1, basis2, ppMethod, ppParams);
outmaps = inmaps;
outmaps.ppMethod = ppMethod;
outmaps.ppParams = ppParams;
if strcmp(ppMethod,'closestRotation')    
    for i=1:length(inmaps.maps)
        for j=1:length(ppParams)
            outmaps.ppMaps{i,j} = closestRotation(inmaps.maps{i}, ppParams{j});
        end        
    end
elseif strcmp(ppMethod,'ICP')    
    for i=1:length(inmaps.maps)
        % 10 iterations
        outmaps.ppMaps{i} = icpFmap(inmaps.maps{i}, basis1, basis2, ppParams{1});
    end
else
    error(['Unrecognized pp method: ' ppMethod]);
end
    
