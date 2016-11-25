function [C,resid,resid_d,resid_o] = leastSquaresMapFast(F,G,S,R,alpha,colMajor)
% If only descriptor constraints:
% Computes a map matrix by solving:
%   min ||CF-G||_2^2
% When more than one such solution exists, return the one with least norm.

if nargin < 5
    alpha = 1;
end
if nargin < 6
    colMajor = 0;
end

% Remember: 
% CF = G <--> F'C' = G', so we solve for C as ...
if isempty(S)
    C = (F'\G')';
    fprintf('Only descriptor constraints used.\n');
% If we have operator constraints, need to flatten constraint matrices
% Solving min ||CF-G||_2^2 + sum_i lambda_i*||R_i C - C S_i||_2^2 + 
% The lambdas are computed using the inf norm of the repsective constraint
% matrices
else
    if colMajor
        [Md,bd] = descriptorConstraints(F,G);
    else
        [Md,bd] = descriptorConstraintsRowMajor(F,G);
    end
    Md = Md/norm(Md,inf);
    bd = bd/norm(Md,inf); 

    M = Md; b = bd;
    for i=1:length(S)
        % Diagonal operators, easier construction
        if isDiag(S{i}) && isDiag(R{i})
            [Mo,bo] = operatorConstraintsDiagonal(S{i},R{i}, colMajor);
        else
            if colMajor
                [Mo,bo] = operatorConstraints(S{i},R{i});
            else
                [Mo,bo] = operatorConstraintsRowMajor(S{i},R{i});
            end
        end
        Mo = alpha*Mo/norm(Mo,inf);
        bo = alpha*bo/norm(Mo,inf); 
        M = [M; Mo];
        b = [b; bo];
    end
    
    % Faster than \ for some reason
    [cs,Rs] = qr(M,b);
    v = Rs\cs;
    
    resid = norm(M*v-b);
    resid_d = norm(Md*v - bd);
    resid_o = norm(Mo*v - bo);
    [mr,nr] = size(R{1});
    [ms,ns] = size(S{1});
    mc = nr; nc = ms;
    if mr ~= mc || nc ~= ns
        error('?');
    end
    if colMajor
        C = full(reshape(v,mc,nc));    
    else
        C = full(reshape(v',nc,mc)');    
    end
    
    %close rotation
%     alpha = 0.01;
%     [uu,ss,vv] = svd(C);
%     dd = diag(ss);
%     
%     dd(dd > alpha) = 1;
%     ss = diag(dd);
%     
%     C = uu*ss*vv';
end