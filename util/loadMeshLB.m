function mesh = loadMeshLB(meshname,varargin)
% This method searches for meshname.mat if it exists and loads that. 
% Otherwise it tries to load meshname.off or meshname.obj, computes LB
% matrix/eigenstuff, and caches it for next time in meshname.mat.  It
% also can load the LB matrix and fill in eigenstuff if necessary.
% loadMeshLB returns two structures used for mesh and LB information.
%
% "Mesh" will be a structure storing all mesh data.

mesh = [];

if length(varargin)>0
    numEigs = varargin{1};
else
    numEigs = -1; % take the default
end


% Check if we've already computed the LB basis etc.
if exist([meshname '_LB_cache.mat'],'file')==2
    load([meshname '_LB_cache.mat']);
    
    fprintf('Loading mesh from cache\n');
    % check if there are enough basis vectors
    if size(mesh.eigenvalues,1) >= numEigs
        % truncate if there are too many basis vectors
        if numEigs > 0
            mesh.eigenvalues = mesh.eigenvalues(1:numEigs);
            mesh.laplaceBasis = mesh.laplaceBasis(:,1:numEigs);
        end
    end
else
    % Load the mesh from .mat file
    if exist([meshname '.mat'],'file')==2
        % Loads a struct called "surface" with elements X, Y, Z, TRIV
        s = load([meshname '.mat']);
       % load('triangles.mat');
    
        %mesh.triangles = tl;
        %mesh.vertices = double(s.shape');
      %mesh.texture = double(s.tex');
        % high-def
        mesh.triangles = s.FV.faces;
        mesh.vertices = s.FV.vertices;
        mesh.texture = s.FV.facevertexcdata;
        
        %mesh.triangles = s.surface.TRIV;
        %mesh.vertices = [s.surface.X s.surface.Y s.surface.Z];
    % Load the mesh from .off file
    else
        [X T] = readOff([meshname '.off']);
        mesh.vertices = X;
        mesh.triangles = T;
    end
    mesh.nv = length(mesh.vertices);
    mesh.nf = length(mesh.triangles);

    % Check if cotangent Laplacian exists and compute if necessary
    if exist([meshname '_cot.mat'],'file')==2
        s = load([meshname '_cot.mat']);
        mesh.areaWeights = s.A;
        mesh.cotLaplace = s.W;
    else
        % Compute LB matrix, using *negative* cot laplacian
        [W A] = cotLaplacian(mesh.vertices, mesh.triangles);
        mesh.areaWeights = 2*A;
        mesh.cotLaplace = -2*W;
    end

    rescaledVerts = mesh.vertices/sqrt(sum(mesh.areaWeights)/2);
    dlmwrite(['cache/' meshname '_rescaled_cloud.xyz'],rescaledVerts,' ');
    
    % % Normalize mesh area to sum to 1
    mesh.areaWeights = mesh.areaWeights / sum(mesh.areaWeights);

    % Compute eigenstuff
    n = size(mesh.vertices,1);
    Am = sparse(1:n,1:n,mesh.areaWeights);
    [evecs, evals] = eigs(mesh.cotLaplace, Am, max(numEigs,100), -1e-5);
    evals = diag(evals);

    mesh.laplaceBasis = evecs;
    mesh.eigenvalues = evals;

    % Compute non-weighted Laplace eigenstuff
    [evecs, evals] = eigs(mesh.cotLaplace, max(numEigs,100), -1e-5);
    mesh.confLaplaceBasis = evecs;
    mesh.confLaplaceEvals = diag(evals);
    
    % Cache eigenstuff so you don't have to recompute
    save(['cache/' meshname '_LB_cache.mat'],'mesh');

    % Truncate eigenstuff if you have too much
    if numEigs > 0
        mesh.eigenvalues = mesh.eigenvalues(1:numEigs);
        mesh.laplaceBasis = mesh.laplaceBasis(:,1:numEigs);
        % confLaplace
        mesh.confLaplaceBasis = mesh.confLaplaceBasis(:,1:numEigs);
        mesh.confLaplaceEvals = mesh.confLaplaceEvals(1:numEigs);
        
    end
end

mesh.name = meshname;
mesh.nv = length(mesh.vertices);
mesh.nf = length(mesh.triangles);
mesh.edgeLens = edgeLengthMatrix(mesh);
