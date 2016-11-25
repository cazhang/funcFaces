function [S,Fmaps,F,G] = symmetryOperator(mesh, options)

mesh1 = mesh;
mesh2 = mesh;

%%%%%% Basis
% use LB basis
basis1 = mesh1.laplaceBasis;
basis2 = mesh2.laplaceBasis;

% Compute descriptors
[f g] = collectDescriptors(mesh1, mesh2, options);

% Write descriptors in provided basis
F = descriptorsToCoeffs(f, basis1);
G = descriptorsToCoeffs(g, basis2);

% Compute the maps
Fmaps = computeMap(basis1, F, basis2, G, options);

S = Fmaps{1}.maps{1};