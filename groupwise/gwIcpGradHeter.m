function g = gwIcpGradHeter(X, Nsets, Ndims, NNbasis)
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
                
                if (i~=j && i~=1)
                    
%                     V1 = V(:,:,j);
%                     V2 = V(:,:,i);                   
%                     Vjj = V1;
%                     Vii = V2;
%                     map = U{i} * T{j};                                                                        
%                     dest1 = annquery(V2, map * V1, 1);
%                     Vji = V2(:, dest1);                  
%                     dest2 = annquery(V1, map' * V2, 1); 
%                     Vij = V1(:, dest2);
                    Vjj = NNbasis{j,j};                   
                    Vji = NNbasis{i,j};
                    Vii = NNbasis{i,i};                   
                    Vij = NNbasis{j,i};
                    tmp1 = Vji - U{i}*T{j}*Vjj;
                    tmp2 = Vij - U{j}*T{i}*Vii;
                    g(:, :, i - 1)  = g(:, :, i - 1) +  T{i}*tmp1 * Vjj'*U{j}*T{i} - T{j} * tmp2 * Vii';
                    clear tmp
                    
%                     % V1 --> V2
%                     dest1 = annquery(V2, map * V1, 1);
%                     % V2 --> V1
%                     dest2 = annquery(V1, map' * V2, 1);      
%                     in = [1:length(V1)];
%                     [~, nnidx] = find(in == dest2(dest1(in)));  
%                     V1 = V1(:, nnidx);
%                     V2 = V2(:, nnidx);
                    
                    
                    
                end
            end
        end
        g = 2 * g;
                
    end