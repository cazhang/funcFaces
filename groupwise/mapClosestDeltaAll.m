function [NNBasis, NNInd] = mapClosestDeltaAll(X, Vin, Nsets, step)


T = M2T(X);
nDim = size(Vin, 1);
N = size(Vin, 2);
nPts = length(1:step:N);
NNInd = zeros(nPts, Nsets, Nsets);
NNBasis = zeros(nDim, nPts, Nsets, Nsets);

newInd = [1:step:N];
%newIndex = randperm(N);
%newInd = newIndex(1:nPts);


for i=1:Nsets
    for j=1:Nsets
        % V1->V2
        V1 = Vin(:,newInd,j);
        V2 = Vin(:,newInd,i);
        
        % from j to i
        M = T(:, :, i, j);
        %R = closestRotation(M);
        Vc =  M * V1;
        
        dest = annquery(V2, Vc, 1);
        
        NNInd(:,i,j) = newInd(dest);
        NNBasis(:, :, i, j) = V2(:, dest);
    end
end


end