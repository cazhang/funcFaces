function G = Gnl_Eigen_Matrix(cds, h, Area)
% This function calculates the Laplacian matrix for the generalized eigenvalue 
% problem. L = G * diag(Area).
% 
% Input:
%   cds         - coordinates of every points. It is a Npts*3 matrix, where
%                 Npts is the number of points, each row in the matrix
%                 corresponds to the coordinate of one point.
%   h           - sqrure of the average distance of any point to its k
%                 nearest neighbors.
%   Area        - The Voronoi area weight for each point, it is a Npts*1 
%                 vector.
% Output:
%   G           - matrix for generalized eigenvalue problem where L = G * diag(Area).
%
% Copyright (c) Chuanjiang Luo, Issam Safa and Yusu Wang, 2009  
% The Ohio State University

Np = size(cds, 1);

Nmax = max(1e+7, floor(Np^2*0.004));
ii(Nmax) = 0;
jj(Nmax) = 0;
kk(Nmax) = 0;

count = 0;

flag = 1 : Np;

for k = 1 : Np
    r2 = 36 * h;
    r = 6 * sqrt(h);
    
    nbindex = flag(abs(cds(:, 1)-cds(k, 1))<=r & abs(cds(:, 2)-cds(k, 2))<=r & abs(cds(:, 3) - cds(k,3)) <= r);
    nbpts = cds(nbindex, :);
    
    s = (nbpts(:, 1)-cds(k, 1)).^2 + (nbpts(:, 2) - cds(k, 2)) .^ 2 + (nbpts(:, 3) - cds(k, 3)) .^ 2;
    ridc = nbindex(s < r2);
    rdist = sqrt(s(s < r2));
   
    digindex = find(rdist == 0);
    
    number = length(ridc);
    dist2 = rdist .^ 2;
    G = -exp(- dist2 ./ (4 * h)) / (4 * pi * h ^ 2);
    
    G(digindex) = -sum(G .* Area(ridc)) / Area(k) + G(digindex);
    
    ii(count + 1: count + number) = k;
    jj(count + 1: count + number) = ridc;
    kk(count + 1: count + number) = G;
    count = count + number;
end

ii(count + 1 : end) = [];
jj(count + 1 : end) = [];
kk(count + 1 : end) = [];

G = sparse(ii, jj, kk, Np, Np);