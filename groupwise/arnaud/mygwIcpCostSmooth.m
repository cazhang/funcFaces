function f = mygwIcpCostSmooth(X, Nsets, Basis, LapMat)
T = myM2T(X);
f = 0;
for i = 1: Nsets
    for j = 1: Nsets
        tmp = T(:, :, i, j) * Basis(:, :, j);
        f   = f + trace(tmp * LapMat{j} * tmp.');
    end
end