function [C,resid,resid_d,resid_o] = leastSquaresMap(F,G,S,R,alpha)
% If only descriptor constraints:
% Computes a map matrix by solving:
%   min ||CF-G||_2^2
% When more than one such solution exists, return the one with least norm.

if nargin < 5
    alpha = 1;
end

% Remember: 
% CF = G <--> F'C' = G', so we solve for C as ...
if isempty(S)
    C = (F'\G')';
% If we have operator constraints, need to flatten constraint matrices
% Solving min ||CF-G||_2^2 + sum_i lambda_i*||R_i C - C S_i||_2^2 + 
% The lambdas are computed using the inf norm of the repsective constraint
% matrices
else
    [Md,bd] = descriptorConstraints(F,G);
    Md = Md/norm(Md,inf);
    bd = bd/norm(Md,inf); 

    M = Md; b = bd;
    for i=1:length(S)
        [Mo,bo] = operatorConstraints(S{i},R{i});
        Mo = alpha*Mo/norm(Mo,inf);
        bo = alpha*bo/norm(Mo,inf); 
        M = [M; Mo];
        b = [b; bo];
    end
    v = M\b;
    resid = norm(M*v-b);
    resid_d = norm(Md*v - bd);
    resid_o = norm(Mo*v - bo);
    [mr,nr] = size(R{1});
    [ms,ns] = size(S{1});
    mc = nr; nc = ms;
    if mr ~= mc || nc ~= ns
        error('?');
    end
    C = full(reshape(v,mc,nc));    
end



    