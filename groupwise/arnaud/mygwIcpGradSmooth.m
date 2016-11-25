function g = mygwIcpGradSmooth(X, Nsets, Ndims, Basis, LapMat)
g   = zeros(Ndims, Ndims, Nsets - 1);
T   = myM2T(X);
T   = T(:, :, 1, :);
TB  = zeros(size(Basis));
for i = 1: Nsets
    TB(:, :, i) = T(:, :, i) * Basis(:, :, i);
end
for i = 2: Nsets
    for j = 1: Nsets
        g(: ,: , i - 1) = g(:, :, i - 1) + TB(:, :, j) * LapMat{j} * (TB(:, :, j).' * T(:, :, i));
    end
    g(: ,: , i - 1) = g(:, :, i - 1) + Nsets * TB(:, :, i) * LapMat{i} * Basis(:, :, i).';
end
g = 2 * permute(g, [2 1 3]);