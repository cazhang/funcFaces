function matrix = svdMap(F, G, alpha)
% Computes a map matrix using orthogonal Procrustes.

xx = F * G'; 
matrix = closestRotation(xx,alpha)';