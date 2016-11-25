function cmaps  = computeMapLinear(M,b,m,n)
% Computes a functional map by solving M\b

cmap.method = 'linear';
v = M\b;
C = reshape(v,m,n);
cmap.maps{1} = closestRotation(full(C),1e-4);
cmaps{1} = cmap;