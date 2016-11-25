function hks = heatKernelSignature(mesh,numTimes)
% This method computes the heat kernel signature for each vertex on a list.
% It uses precomputed LB eigenstuff stored in "mesh" and automatically
% chooses the time steps based on mesh geometry.

numVertices = size(mesh.vertices, 1);

tmin = -4*log(10)/mesh.eigenvalues(end);
tmax = -4*log(10)/mesh.eigenvalues(10); % Why 10?
ts = logspace(log10(tmin),log10(tmax),numTimes);

areaMtx = sparse(1:numVertices,1:numVertices,mesh.areaWeights);
D = mesh.laplaceBasis' * (areaMtx * mesh.laplaceBasis.^2);

T = exp(-abs(mesh.eigenvalues)*ts);
F = D*T;

fs = mesh.laplaceBasis * F;    
nf = size(F,2);
scale = sparse(1:nf, 1:nf, 1./max(fs), nf, nf);
fs = fs*scale;    
hks =  mesh.laplaceBasis'*(areaMtx*fs);

hks = mesh.laplaceBasis * hks;