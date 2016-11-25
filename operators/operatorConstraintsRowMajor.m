function [M,b] = operatorConstraintsRowMajor(Sf,Rf)
% Given the functional operators Sf on mesh1 and Rf on mesh2, build the linear
% constraints such that Mc = 0 is equivalent to RC = CS, where 
% c = reshape(C',size(C,2)*size(C,1),1)

[mr,nr] = size(Rf);
[ms,ns] = size(Sf);

% operators should be square
if mr ~= nr || ms ~= ns
    error('operators should be square');
end

mc = nr;
nc = ms;

I = zeros(mc*nc*(mc+nc),1);
J = I; S = I;
last = 1;
for i=1:mc
    for j=1:nc
        % Equation for the [i,j]th entry of the matrix = row ii in M
        ii = sub2ind_r([mc,nc],i,j);
        % R part 
        % j-th column of C
        Jr = sub2ind_r([mc,nc],(1:mc),repmat(j,1,mc));
        % i-th row of R 
        Sr = Rf(i,:);

        % S part 
        % i-th row of C
        Js = sub2ind_r([mc,nc],repmat(i,1,nc),(1:nc));
        % j-th column of S
        Ss = -Sf(:,j)';
        
        % Indices for current element
        Jc = [Jr,Js];
        Sc = [Sr,Ss];
        Ic = repmat(ii,size(Jc));
        
        I(last:last+length(Ic)-1) = Ic;
        J(last:last+length(Jc)-1) = Jc;
        S(last:last+length(Sc)-1) = Sc;
        
        last = last + length(Ic);
    end
end

M = sparse(I,J,S,mc*nc,mc*nc);
b = sparse(size(M,2),1);