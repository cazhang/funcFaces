
function g = gwIcpGradConsisSmooth(X, Nsets, Ndims, NNbasis, NNIndGPA, Basis, LapMat, w1, w2, w3)
g = zeros(Ndims, Ndims, Nsets-1);
T = cell(Nsets, 1);
U = cell(Nsets, 1);
TT = M2T(X);
for i=1:Nsets
    T{i} = TT(:,:,1,i);
    U{i} = TT(:,:,i,1);
end
nDim = size(Basis, 1);
% add the gradient for temp1
for i = 1:Nsets
    for j = 1:Nsets
        
        if (i~=j && i~=1)
            map = U{i}*T{j};
            % ICP
            VVjj = NNbasis(:,:,j,j);
            VVji = NNbasis(:,:,i,j);
            VVii = NNbasis(:,:,i,i);
            VVij = NNbasis(:,:,j,i);
            tmp1 = VVji - U{i}*T{j}*VVjj;
            tmp2 = VVij - U{j}*T{i}*VVii;
            g(:, :, i - 1)  = g(:, :, i - 1) +  (w2^2)*(T{i}*tmp1 * VVjj'*U{j}*T{i} - T{j} * tmp2 * VVii');
            clear tmp1 tmp2
            
            % landmark
            Vjj = Basis(:,NNIndGPA(:,j),j);
            Vji = Basis(:,NNIndGPA(:,i),i);
            Vii = Basis(:,NNIndGPA(:,i),i);
            Vij = Basis(:,NNIndGPA(:,j),j);

            tmp1 = Vji - U{i}*T{j}*Vjj;
            tmp2 = Vij - U{j}*T{i}*Vii;
            g(:, :, i - 1)  = g(:, :, i - 1) +  (w1^2)*(T{i}*tmp1 * Vjj'*U{j}*T{i} - T{j} * tmp2 * Vii');
            clear tmp1 tmp2
            
            % smoothness
            
            %tmp1 = U{i}' * (map*Basis(:,:,j)) * LapMat{j} * (map*Basis(:,:,j))';
            %tmp2 = U{i}' * (map * Basis(:,:,j)) * LapMat{j} * Basis(:,:,j)';
            %g(:,:,i-1) = g(:,:,i-1) + (w3)*(tmp2-tmp1);
            %clear tmp1 tmp2
%             for k = 1:nDim
%                 P1 = T{j}*Basis(:,:,j)*LapMat{j}*Basis(:,:,j)'*T{j}';
%                 P2 = Basis(:,:,j)*LapMat{j}*Basis(:,:,j)';
%                 
%                 tmp1 = (P1+P1')*T{i};
%                 tmp2 = (P2+P2')*T{j}'*T{i}*T{i};
%                 g(:,:,i-1) = g(:,:,i-1) + (w3^2)*(tmp1+tmp2);
%                 clear tmp1 tmp2
%             end         
        end
    end
end
T   = TT(:, :, 1, :);
TB  = zeros(size(Basis));
for i = 1: Nsets
    TB(:, :, i) = T(:, :, i) * Basis(:, :, i);
end
for i = 2: Nsets
    for j = 1: Nsets
        g(: ,: , i - 1) = g(:, :, i - 1) + TB(:, :, j) * LapMat{j} * (TB(:, :, j).' * T(:, :, i));
    end
    g(: ,: , i - 1) = g(:, :, i - 1) + Nsets * TB(:, :, i) * LapMat{i} * Basis(:, :, i).';
end
g = g + permute(g, [2 1 3]);
            
g = 2 * g;

end