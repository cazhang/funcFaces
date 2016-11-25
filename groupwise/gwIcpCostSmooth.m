function f = gwIcpCostSmooth(X, Nsets, Basis, LapMat)
        f = 0;
        
        % mapping from i to 1
        T = cell(Nsets, 1);
        % mapping from 1 to i
        U = cell(Nsets, 1);
        % generate all Nsets^2 mapping from (Nsets-1) maps
        TT = M2T(X);
        for i=1:Nsets
            T{i} = TT(:,:,1,i);
            U{i} = TT(:,:,i,1);
        end
      
        
        for i = 1:Nsets
            for j = 1:Nsets
                    if (i~=j)
                    % mapping from j to i
                    %map = U{i} * T{j};
                 
                    % Smoothness term
                    %I = ones(nDim, 1);
                    %t1 = map * Basis(:,:,j) * LapMat{j} * Basis(:,:,j)' * map';
                    %temp1 = trace(t1); % ignore the 2 times, to prevent multiply 2 in grad.
                    %f = f + w*temp1;
                    %map = map';
                    %t2 = map * Basis(:,:,i) * LapMat{i} * Basis(:,:,i)' * map';
                    %temp2 = trace(t2); % ignore the 2 times, to prevent multiply 2 in grad.
                    tmp = (U{i} * T{j}) * Basis(:,:,j);

                    f = f + trace(tmp * LapMat{j} * tmp');
                                                    
                    
                    
                    clear tmp
                    end
                
            end
        end
        
    end