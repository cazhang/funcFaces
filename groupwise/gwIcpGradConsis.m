function g = gwIcpGradConsis(X, Nsets, Ndims, NNbasis, NNIndGPA, Basis, w)
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
        
        % for the w_star
        % mapping from j to i
        map = U{i} * T{j};
        % Landmark term
        V1 = Basis(:,NNIndGPA(:,j),j);
        V2 = Basis(:,NNIndGPA(:,i),i);
        
        % ICP term
        VV1 = NNbasis(:,:,j,j);
        VV2 = NNbasis(:,:,i,j);
        
        temp1 = V2 - map * V1; % landmark
        temp2 = VV2 - map * VV1; % nn term
        
        norm_nn = norm(temp2(:),2)^2;
        norm_lm = norm(temp1(:),2)^2;
        
        ratio = norm_lm / norm_nn;
        
        w_star = w/ratio;
                    
        
        if (i~=j && i~=1)
            map = U{i}*T{j};
            % ICP
            VVjj = NNbasis(:,:,j,j);
            VVji = NNbasis(:,:,i,j);
            VVii = NNbasis(:,:,i,i);
            VVij = NNbasis(:,:,j,i);
            tmp1 = VVji - U{i}*T{j}*VVjj;
            tmp2 = VVij - U{j}*T{i}*VVii;
            g(:, :, i - 1)  = g(:, :, i - 1) +  (T{i}*tmp1 * VVjj'*U{j}*T{i} - T{j} * tmp2 * VVii');
            clear tmp1 tmp2
            
            % landmark
            Vjj = Basis(:,NNIndGPA(:,j),j);
            Vji = Basis(:,NNIndGPA(:,i),i);
            Vii = Basis(:,NNIndGPA(:,i),i);
            Vij = Basis(:,NNIndGPA(:,j),j);

            tmp1 = Vji - U{i}*T{j}*Vjj;
            tmp2 = Vij - U{j}*T{i}*Vii;
            g(:, :, i - 1)  = g(:, :, i - 1) +  (w^2) * (T{i}*tmp1 * Vjj'*U{j}*T{i} - T{j} * tmp2 * Vii');
            clear tmp1 tmp2
                    
        end
    end
end
g = 2 * g;

end