function [M,b] = descriptorConstraintsRowMajor(F,G)
% Build the linear constraints such that Mc = 0 is equivalent to CF = G, where 
% c = reshape(C',size(C,2)*size(C,1),1)

[mf,nf] = size(F);
[mg,ng] = size(G);

% Num columns should be the same
if nf ~= ng
    error('num columns should be the same');
end

mc = mg;
nc = mf;

% Sorry about this... It's much faster this way 
if mc == nc
    % Saves the sort but will probably not work for non square C
    I = reshape((repmat([1:nf],mc,nc) + repmat(((0:mc-1)*nf)',1,mc*nf))',mc*nc*nf,1);
    J = reshape(repmat(1:nc*mc,nf,1),nc*mc*nf,1);
    S = repmat(reshape(F',nf*mf,1),mc,1);
else
    I = reshape(repmat((1:mc*nf),nc,1),mc*nf*nc,1);
    J = reshape((repmat(1:nc,mc,nf) + repmat((0:mc:(mc-1)*mc)',1,nf*nc))',1,mc*nf*nc)';
    S = reshape(repmat(reshape(F,1,size(F,2)*size(F,1)),mc,1)',1,mc*nf*nc)';

    % Faster, since matlab sparse storage is column major
    [~,is] = sort(J);
    J=J(is);
    I=I(is);
    S=S(is);
end

M = sparse(I,J,S,mc*nf,mc*nc);
b = reshape(G', size(G,2)*size(G,1),1);

