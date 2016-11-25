function matrix = groundTruthMap(basis1, basis2)
% Computes a ground truth map matrix assuming the two meshes are in
% vertex-to-vertex correspondence.

% We assume that the meshes have the same vertices, so functions get
% carried over directly.  Thus, we just solve a least-squares problem:
matrix = basis2\basis1;