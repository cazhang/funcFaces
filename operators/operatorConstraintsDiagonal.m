function [M,b] = operatorConstraintsDiagonal(Sf,Rf,colMajor)
% Given the functional operators Sf on mesh1 and Rf on mesh2, build the linear
% constraints such that Mc = 0 is equivalent to RC = CS, where 
% colMajor = 1 ==> c = reshape(C,size(C,2)*size(C,1),1)
% colMajor = 0 ==> c = reshape(C',size(C,2)*size(C,1),1)
% Assumes Sf and Rf are diagonal

[mr,nr] = size(Rf);
[ms,ns] = size(Sf);

% operators should be square
if mr ~= nr || ms ~= ns
    error('operators should be square');
end

if ~isDiag(Sf) || ~isDiag(Rf)
    error('operators should be diagonal');
end

[I,J] = meshgrid(1:ms,1:mr);
dS = diag(Sf); dR = diag(Rf);
S = dS(I) - dR(J);

nM = size(S,1)*size(S,2);
if colMajor
    Sflat = reshape(S, nM, 1);
else
    Sflat = reshape(S', nM,1);
end

M = spdiags(Sflat,0,nM,nM);
b = sparse(size(M,2),1);