function [NNBasis] = mapClosestDeltaAllHeter(X, Vin, Nsets, step)


T = M2T(X);
nDim = size(Vin{1,1}, 1);
temp = [];
for i=1:Nsets
    temp = [temp size(Vin{1,i}, 2)];
end
N = min(temp);
nPts = length(1:step:N);
NNBasis = cell(Nsets, Nsets);

newInd = [1:step:N];
%newIndex = randperm(N);
%newInd = newIndex(1:nPts);


for i=1:Nsets
    for j=1:Nsets
        % V1->V2
        V1 = Vin{1,j}(:,newInd);
        V2 = Vin{1,i}(:,newInd);
        
        % from j to i
        M = T(:, :, i, j);
        %R = closestRotation(M);
        Vc =  M * V1;
        
        dest = annquery(V2, Vc, 1);
        NNBasis{i, j} = V2(:, dest);
    end
end


end