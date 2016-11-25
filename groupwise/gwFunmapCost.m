function f = gwFunmapCost(X, P, Q, Nsets)
        f = 0;
        T = cell(Nsets, 1);
        U = cell(Nsets, 1);
        TT = M2T(X);
        for i=1:Nsets
            T{i} = TT(:,:,1,i);
            U{i} = TT(:,:,i,1);
        end
        
        for i = 1: Nsets
            for j = 1: Nsets
                if i ~= j
                    map = U{i} * T{j};
                    % descriptor constraints
                    temp1 = P(:, :, i) - map * P(:, :, j);
                    %temp1 = temp1 / norm(temp1, inf);             
                    % operator constraints
                    temp2 = Q(:, :, j) * map - map * Q(:, :, i);
                    %temp2 = temp2 / norm(temp2, inf);
                    
                    % descriptor only
                    temp = [temp1(:); temp2(:)];
                    
                    f = f + norm(temp, 2)^2;
                    
                    clear temp1 temp2 temp
                end
            end
        end
    end