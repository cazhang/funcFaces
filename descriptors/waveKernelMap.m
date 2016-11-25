function wkm = waveKernelMap(mesh,numTimes,vertex)
% This method computes the wave kernel map for each vertex on a list.  From
% the given vertex parameter.
% It uses precomputed LB eigenstuff stored in "mesh" and automatically
% chooses the time steps based on mesh geometry.

numEigenfunctions = size(mesh.eigenvalues,1);

absoluteEigenvalues = abs(mesh.eigenvalues);
emin = log(absoluteEigenvalues(2));
emax = log(absoluteEigenvalues(end));
s = 7*(emax-emin) / numTimes; % Why 7?
emin = emin + 2*s;
emax = emax - 2*s;
es = linspace(emin,emax,numTimes);

T = exp(-(repmat(log(absoluteEigenvalues),1,numTimes) - repmat(es,numEigenfunctions,1)).^2/(2*s^2));
wkm = T.*repmat(mesh.laplaceBasis(vertex,:)', 1, numTimes);
wkm = mesh.laplaceBasis * wkm;