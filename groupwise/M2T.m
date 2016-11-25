function T = M2T(M)
[N, N, n] = size(M);
T   = zeros(N, N, n+1, n+1);
% Prepare Transform matrix
for i = 1: n+1
    if i == 1
        T(:,:,1,i) = eye(N);
    else
        T(:,:,1,i) = M(:,:,i-1);
    end
    
end

for i = 2: n+1
    T(:, :, i, i)   = eye(N);
    %T(:, :, 1, i)   = T(:, :, i, 1)^-1;
    T(:, :, i, 1) = T(:, :, 1, i).';
    for j = 2: i - 1
        T(:, :, i, j) = T(:, :, i, 1) * T(:, :, 1, j);
        T(:, :, j, i) = T(:, :, j, 1) * T(:, :, 1, i);
    end
    
end

end