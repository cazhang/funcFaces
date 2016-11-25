function f = gwIcpCostConsis(X, Nsets, NNbasis, NNIndGPA, Basis, w)
        f = 0;
        
        % mapping from i to 1
        T = cell(Nsets, 1);
        % mapping from 1 to i
        U = cell(Nsets, 1);
        TT = M2T(X);
        for i=1:Nsets
            T{i} = TT(:,:,1,i);
            U{i} = TT(:,:,i,1);
        end
        nDim = size(Basis, 1);
        
        for i = 1:Nsets
            for j = 1:Nsets
                if i ~= j
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
                    
                    %disp(w_star);            
                    
                    temp1 = w * temp1;
                    % descriptor only
                    %temp = [temp1(:); temp2(:)];
                    
                    f = f + norm(temp1(:), 2)^2 + norm(temp2(:), 2)^2;
                    
                    clear temp1 temp2 temp
                end
            end
        end
    end