% This script calculates the eigen-decomposition of the Voronoi-Laplace
% matrix by the generalized eigenvalue problem. The Laplace matrix L = G * D;
%
% Copyright (c) Chuanjiang Luo, Issam Safa and Yusu Wang, 2009  
% The Ohio State University

addpath([cd '\functions']);
addpath([cd '\data']);
addpath([cd '\GLtree3DMex']);
%matlabpool

filename = 'sphere_2562.pcd';       % set to the file of the model. 
cds = load(filename);
N = 100;                            % set the number of eigenvalue you want to calculate here
tsn = 5;  %was 15                                     % set the number of points to approximate the tangent plane here
tic
Npts = size(cds, 1);

ptrtree = BuildGLTree(cds);
[kidc,kdist] = KNNSearch(cds, cds, ptrtree, tsn);
DeleteGLTree(ptrtree);
kidc = int32(kidc);

% end-10:end-1 -- old subscripts on kdist
h = (sum(sum(kdist)) / 10 / Npts) ^ 2;         % Note if you want to compare the eigenvalue of two shapes, you need to set h the same for both shapes.

AW = VoronoiArea(cds, kidc, sqrt(h));

G = Gnl_Eigen_Matrix(cds, h, AW);

invD = sparse(1 : Npts, 1 : Npts, 1 ./ AW);
[eigenvector, eigenvalue] = eigs(G, invD, N, -1e-5);
%clear G
eigenvalue = diag(eigenvalue);

if eigenvalue(1) > eigenvalue(end)
    eigenvalue = eigenvalue(N : -1 : 1);
    eigenvector = eigenvector(:, N : -1 : 1);
end

for i = 1 : N
    eigenvector(:, i) = eigenvector(:, i) ./AW;
end

toc

%matlabpool close