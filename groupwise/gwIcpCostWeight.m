function f = gwIcpCostWeight(X, Nsets, NNbasisGPA, NNIndGPA, NNbasis, w)
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
        
        for i = 1:Nsets
            for j = 1:Nsets
                if i ~= j
                    % mapping from j to i
                    map = U{i} * T{j};
                    % Mutual NN
                    V1 = NNbasisGPA(:,NNIndGPA{j,j},j);
                    V2 = NNbasisGPA(:,NNIndGPA{i,j},i);
                    
                    % Single NN
                    VV1 = NNbasis(:,:,j,j);
                    VV2 = NNbasis(:,:,i,j);
                    
                    %V1 = V(:,:,j);
                    %V2 = V(:,:,i);                   
                    %Vc = map * V1;
                    %dest = annquery(V2, Vc, 1);                   
                    %V2 = V2(:, dest);
                                                     
%                     % mapped basis of j
%                     % ann search
%                     % V1 --> V2
%                     dest1 = annquery(V2, map * V1, 1);
%                     V2 = V2(:, dest1);
%                     % V2 --> V1
%                     dest2 = annquery(V1, map' * V2, 1);
%                     
%                     % for i-th point in V1, the est. corres in V2 is
%                     % dest1(i). and the corres of this point in V1 is
%                     % dest2(dest1(i))
%                     in = [1:length(V1)];
%                     [~, nnidx] = find(in == dest2(dest1(in)));                   
%                     %nnidx = annquery(V2, Vc, 1);
%                     V1 = V1(:, nnidx);
%                     V2 = V2(:, nnidx);
                    
                    % distance in functional space
                    % is this map necessary??
                    temp1 = V2 - map * V1;
                    temp2 = VV2 - map * VV1;
                    temp1 = temp1 * w;
                    
                    
                    
                    % descriptor only
                    temp = [temp1(:); temp2(:)];
                    
                    f = f + norm(temp, 2)^2;
                    
                    clear temp1 temp2 temp
                end
            end
        end
    end