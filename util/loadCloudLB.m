function cloud = loadCloudLB(cloudname,numEigenfunctions)

spacing = 1;

fprintf(1,'Loading cloud...');
cloud.vertices = load([cloudname '.xyz']);
cloud.vertices = cloud.vertices(1:spacing:size(cloud.vertices,1),:);
cloud.nv = size(cloud.vertices,1);
cloud.nf = 0;

tangentPoints = 10; % default is 15
fprintf(1,'Number of points to approximate tangent plane: %d\n', tangentPoints);

fprintf(1,'KNN search...\n');
ptrtree = BuildGLTree(cloud.vertices);
[kidc,kdist] = KNNSearch(cloud.vertices, cloud.vertices, ptrtree, tangentPoints);
DeleteGLTree(ptrtree);

kidc = int32(kidc);

% This is an approximation of the average edge length; not sure what the
% factor of 10 is doing there -- from example code.
h = (sum(sum(kdist)) / 10 / cloud.nv) ^ 2;
fprintf(1,'Average edge length approx. = %g\n',h);

%matlabpool
fprintf(1,'Computing areas...\n');
AW = VoronoiArea(cloud.vertices, kidc, sqrt(h));

fprintf(1,'Total area: %g\n', sum(AW));
%AW = AW / sum(AW); % Normalize so that shape has area 1 (does this work?)

% Throw out vertices with zero area weight (why does this happen?)
nonzeros = find(AW);
AW = AW(nonzeros);
cloud.vertices = cloud.vertices(nonzeros,:);
cloud.nv = size(cloud.vertices,1);

% Laplace = G * D
fprintf(1,'Computing Laplacian matrix...\n');
G = Gnl_Eigen_Matrix(cloud.vertices, h, AW);
%matlabpool close

cloud.cotLaplace = G;
cloud.areaWeights = AW;

numEigs = numEigenfunctions;
Npts = cloud.nv;

fprintf(1,'Computing %d eigenvalues/eigenvectors...\n', numEigs);
invD = sparse(1 : Npts, 1 : Npts, 1 ./ AW);
[eigenvector, eigenvalue] = eigs(G, invD, numEigs, -1e-5);
eigenvalue = diag(eigenvalue);

if eigenvalue(1) > eigenvalue(end)
    eigenvalue = eigenvalue(numEigs : -1 : 1);
    eigenvector = eigenvector(:, numEigs : -1 : 1);
end

for i = 1 : numEigs
    eigenvector(:, i) = eigenvector(:, i) ./AW;
end

cloud.laplaceBasis = eigenvector;
cloud.eigenvalues = -eigenvalue;
cloud.cotLaplace = -cloud.cotLaplace;
cloud.name = cloudname;
cloud.triangles = [];

%scatter3(cloud.vertices(:,1),cloud.vertices(:,2),cloud.vertices(:,3),3,cloud.laplaceBasis(:,2),'filled');axis equal;axis off;colorbar;