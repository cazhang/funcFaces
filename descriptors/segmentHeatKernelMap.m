function hkm = segmentHeatKernelMap(mesh,numTimes,segment)
% This method computes the heat kernel map for the segment indicator.  From
% the given vertex parameter.
% It uses precomputed LB eigenstuff stored in "mesh" and automatically
% chooses the time steps based on mesh geometry.

tmin = -4*log(10)/mesh.eigenvalues(end);
tmax = -4*log(10)/mesh.eigenvalues(10); % Why 10?
ts = logspace(log10(tmin),log10(tmax),numTimes);

T = exp(-abs(mesh.eigenvalues)*ts);
F = T.*repmat(mesh.laplaceBasis'*segment, 1, length(ts));

fs = mesh.laplaceBasis*F;
numVertices = size(mesh.vertices, 1);
areaMtx = sparse(1:numVertices,1:numVertices,mesh.areaWeights);
hkm =  mesh.laplaceBasis'*(areaMtx*fs);
hkm = mesh.laplaceBasis*hkm;