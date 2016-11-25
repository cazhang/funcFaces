function [NNBasis, NNInd] = mapClosestDeltaDouble(X, Vin, Nsets, step)


T = M2T(X);
nDim = size(Vin, 1);
N = size(Vin, 2);
nPts = length(1:step:N);
NNInd = cell( Nsets, Nsets);
NNBasis = zeros(nDim, nPts, Nsets, Nsets);

for i=1:Nsets
    NNBasis(:,:,i,i) = Vin(:,1:step:end,i);
end

newInd = [1:step:N];


for i=1:Nsets
    for j=1:Nsets
        if i >= j
            % V1->V2
            V1 = Vin(:,1:step:end,j);
            V2 = Vin(:,1:step:end,i);
            
            % from j to i
            M = T(:, :, i, j);
            %R = closestRotation(M);
            Vc =  M * V1;
            dest1 = annquery(V2, Vc, 1);
            % from i to j
            M = T(:, :, j, i);
            %R = closestRotation(M);
            Vc =  M * V2;
            dest2 = annquery(V1, Vc, 1);
            
            % V1 -> V2 (j -> i)
            in = [1:nPts];
            [~, localidx] = find(in == dest2(dest1(in)));
            nnidx = newInd(localidx);
            NNInd{i,j}.src = nnidx;
            nnidx = dest1(localidx);
            NNInd{i,j}.des = newInd(nnidx);
            % V2 -> V1 (i -> j)
            [~, localidx] = find(in == dest1(dest2(in)));
            nnidx = newInd(localidx);
            NNInd{j,i}.src = nnidx;
            nnidx = dest2(localidx);
            NNInd{j,i}.des = newInd(nnidx);
            %nSel = length(nnidx);
            %NNBasis(:, nnidx, i, j) = V2(:, dest1(nnidx));
            %NNBasis(:, nnidx, j, i) = V1(:, dest2(nnidx));
        end
    end
end


end