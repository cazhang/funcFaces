function [M,b] = descriptorConstraints(F,G)
% Build the linear constraints such that Mc = 0 is equivalent to CF = G, where 
% c = reshape(C,size(C,1)*size(C,2),1)

[mf,nf] = size(F);
[mg,ng] = size(G);

% Num columns should be the same
if nf ~= ng
    error('num columns should be the same');
end

mc = mg;
nc = mf;

I = zeros(mc*nf*nc,1);
J = I; S = I;
last = 1;
for i=1:mc
    for j=1:nf
        % Equation for the [i,j]th entry of the matrix = row ii in M
        ii = sub2ind([mc,nf],i,j);
        % F part 
        % i-th row of C
        Js = sub2ind([mc,nc],repmat(i,1,nc),(1:nc));
        % j-th column of F
        Ss = F(:,j)';
       
        % Indices for current element
        Jc = Js;
        Sc = Ss;
        Ic = repmat(ii,size(Jc));
        
        I(last:last+length(Ic)-1) = Ic;
        J(last:last+length(Jc)-1) = Jc;
        S(last:last+length(Sc)-1) = Sc;
        
        last = last + length(Ic);
    end
end

M = sparse(I,J,S,mc*nf,mc*nc);
b = reshape(G, size(G,1)*size(G,2),1);