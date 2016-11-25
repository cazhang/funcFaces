function g = gwIcpGradSmooth(X, Nsets, Ndims, Basis, LapMat)
g = zeros(Ndims, Ndims, Nsets-1);
T = cell(Nsets, 1);
U = cell(Nsets, 1);
TT = M2T(X);
for i=1:Nsets
    T{i} = TT(:,:,1,i);
    U{i} = TT(:,:,i,1);
end
% add the gradient for temp1
for i = 1:Nsets
    for j = 1:Nsets   
        if (i~=1 && i~=j)
            map = U{i}*T{j};                             
            % Arnaud's
            %tmp1 = (T{i}'*T{j}*Basis(:,:,j)) * LapMat{j} * (T{i}'*T{j}*Basis(:,:,j))'*T{i};
            %tmp2 = (T{j}'*T{i}*Basis(:,:,i)) * LapMat{i} * (T{j}'*Basis(:,:,i))';
       
            % Chao's (T^{-1} = T^{T} simplification)
            tmp1 = T{j}*Basis(:,:,j)*LapMat{j}*Basis(:,:,j)'*T{j}'*T{i};
            tmp2 = T{j}*T{j}'*T{i}*Basis(:,:,i)*LapMat{i}*Basis(:,:,i)';
            %map = map';
            %tmp3 = U{j}' * (map*Basis(:,:,i)) * LapMat{i} * (map*Basis(:,:,i))';
            %tmp4 = U{i}' * (map'*Basis(:,:,j)) * LapMat{j} * Basis(:,:,j)';
            g(:,:,i-1) = g(:,:,i-1) + (tmp2 + tmp1);
            clear tmp1 tmp2
                   
        end
    end
end
g = 2 * g;

end