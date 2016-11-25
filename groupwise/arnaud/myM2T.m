function T = myM2T(M)
[N, N, n]           = size(M);
T                   = zeros(N, N, n + 1, n + 1);
T(:, :, 1, 1)       = eye(N);
T(:, :, 1, 2: end)  = M;
T(:, :, 2: end, 1)  = permute(M, [2 1 3]);
for i = 2: n + 1
    T(:, :, i, i) = eye(N);
    for j = 2: i - 1
        T(:, :, i, j) = T(:, :, i, 1) * T(:, :, 1, j);
        T(:, :, j, i) = T(:, :, i, j).';
    end
end