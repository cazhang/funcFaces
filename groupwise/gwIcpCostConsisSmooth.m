function f = gwIcpCostConsisSmooth(X, Nsets, NNbasis, NNIndGPA, Basis, LapMat, w1, w2, w3)
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
                    
                    % Smoothness term
                    %I = ones(nDim, 1);
                    t3 = map * Basis(:,:,j) * LapMat{j} * Basis(:,:,j)' * map';
                    temp3 = trace(t3); % ignore the 2 times, to prevent multiply 2 in grad.
           
          
                    temp1 = V2 - map * V1;
                    temp2 = VV2 - map * VV1;
                    temp1 = temp1 * w1;
                    temp2 = temp2 * w2;
                    temp3 = temp3 * w3;
                    
                    
                    
                    % descriptor only
                    temp = [temp1(:); temp2(:)];
                    
                    f = f + norm(temp, 2)^2 + temp3;
                    
                    clear temp1 temp2 temp3 temp
                end
            end
        end
    end