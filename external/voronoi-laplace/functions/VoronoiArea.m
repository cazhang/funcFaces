function A = VoronoiArea(cds, kidc, ave)
% A = VoronoiArea(nk, cds)
% 
% This function calculates the Voronoi area weight for point clouds data.
% 
% Input:
%   cds         - coordinates of every points. It is a Npts*3 matrix, where
%                 Npts is the number of points, each row in the matrix
%                 corresponds to the coordinate of one point.
%   kidc        - index of 200 nearest neighbour for each points.
%   ave         - estimated average edge length.
% Output:
%   A           - The Voronoi area weight for each point, it is a Npts*1 
%                 vector.
%
% Copyright (c) Chuanjiang Luo, Issam Safa and Yusu Wang, 2009  
% The Ohio State University

Np = size(cds, 1);

A(Np) = 0;
r2 = 36 * ave ^ 2;
r = 6 * ave;
flag = 1 : Np;

parfor k = 1 : Np
    ngbindex = kidc(k, :);
    ngbpts = cds(ngbindex, :);            
    coefs = princomp(ngbpts);
    
    nbindex = flag(abs(cds(:, 1)-cds(k, 1))<= r & abs(cds(:, 2)-cds(k, 2))<= r & abs(cds(:, 3) - cds(k, 3)) <= r);
    nbpts = cds(nbindex, :);
    
    s = (nbpts(:, 1)-cds(k, 1)).^2 + (nbpts(:, 2) - cds(k, 2)) .^ 2 + (nbpts(:, 3) - cds(k, 3)) .^ 2;
    ngbindex = nbindex(s < r2);
    rdist = s(s < r2);
    
    ngbpts = cds(ngbindex, :);
   
    x = ngbpts * coefs(:, 1);
    y = ngbpts * coefs(:, 2);
    scores = [x, y];
    [V, C] = voronoin(scores);
    XX = V(C{rdist == 0}, :);
    if any(isinf(XX(:)))
        [c1, c2] = find(isinf(XX)); %#ok<NASGU>
        XX(c1(1), :) = [];
    end
    
    if size(XX, 1) >= 3
        [K, A(k)] = convhulln(XX); 
    else
        A(k) = 0;
    end
end
A = A';