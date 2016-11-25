function g = gwFunmapGrad(X, P, Q, Nsets, Ndims)
g = zeros(Ndims, Ndims, Nsets-1);
T = cell(Nsets, 1);
U = cell(Nsets, 1);
TT = M2T(X);
for i=1:Nsets
    T{i} = TT(:,:,1,i);
    U{i} = TT(:,:,i,1);
end

% add the gradient for temp1
for i=1:Nsets
    for j = 1:Nsets
        
        if (i~=j && i~=1)
            map = U{i} * T{j};
            % descriptor constraints
            temp1 = P(:, :, i) - map * P(:, :, j);
            %norm1 = 1 / norm(temp1, inf);
            % operator constraints
            temp2 = Q(:, :, j) * map - map * Q(:, :, i);
            %norm2 = 1 / norm(temp2, inf);
            norm1 = 1;
            norm2 = 1;
            tmp             = (U{i} * T{j}) * P(:, :, j);
            g(:, :, i - 1)  = g(:, :, i - 1) + norm1 * (U{i}.' * (P(:, :, i) - tmp) * tmp.' - U{j}.' * (P(:, :, j) - (U{j} * T{i}) * P(:, :, i)) * P(:, :, i).');
            clear tmp
            % add the gradient for temp2
            tmp = (Q(:,:,i) * U{j} * T{i}) - (U{j} * T{i} * Q(:,:,j));
            tmp2 = (Q(:,:,j) * U{i} * T{j}) - (U{i} * T{j} * Q(:,:,i));
            g(:, :, i-1) = g(:, :, i-1) + norm2 * ( (U{j}.' * Q(:,:,i).' * tmp ) - (U{j}.' * tmp * Q(:,:,j).') );
            g(:, :, i-1) = g(:, :, i-1) + norm2 * ( (U{i}.' * tmp2 * Q(:,:,i).' * T{j}.' * U{i}.')-(U{i}.' * Q(:,:,j).' * tmp2 * T{j}.' * U{i}.') );
            clear tmp tmp2
        end
    end
end
g = 2 * g;

end