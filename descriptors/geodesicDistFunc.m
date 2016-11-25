function dist = geodesicDistFunc( mesh, vertex )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
X = mesh.vertices';
F = mesh.triangles';
[dist,S1,Q1] = perform_fast_marching_mesh(X, F, vertex);

end

