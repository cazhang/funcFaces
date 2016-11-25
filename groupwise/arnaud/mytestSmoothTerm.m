% test smooth term
clear all
close all
clc;
% Add path to Manopt
addpath(genpath('./manopt'));

% Add path to ConstructW
addpath('./utils');

% generate X and V
Nsets = 5; % number of objects
Ndims = 30;% number of dimenstion
Npts = 100;% number of points
V = rand(Ndims, Npts, Nsets);
X = rand(Ndims, Ndims, Nsets-1);
% generate W and L
L = [];
w_options = [];
w_options.NeighborMode = 'KNN';
w_options.k = 5; % neet to be tuned
w_options.WeightMode = 'HeatKernel';
w_options.t = 1;

for i=1:Nsets
    W{i} = constructW(V(:,:,i)', w_options);
    L{i} = diag(sum(W{i})) - W{i};
end
    
manifold = stiefelfactory(Ndims, Ndims, Nsets-1);
gwIcpProblem.M = manifold;

gwIcpProblem.cost = @(X) mygwIcpCostSmooth(X, Nsets, V, L);
gwIcpProblem.grad = @(X) gwIcpProblem.M.egrad2rgrad(X, mygwIcpGradSmooth(X, Nsets, Ndims, V, L));
            
checkgradient(gwIcpProblem);