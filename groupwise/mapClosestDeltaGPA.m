function [NNBasis, NNInd] = mapClosestDeltaGPA(X, Vin, Nsets, step)


T = M2T(X);
nDim = size(Vin, 1);
N = size(Vin, 2);
nPts = length(1:step:N);
%NNInd = cell( Nsets, Nsets);
modelSpan = zeros(Nsets+1,1);
newInd = [1:step:N];

numOfPoints = 0;
NNBasis = zeros(nDim, nPts, Nsets, Nsets);
NNInd = cell(Nsets, Nsets);

for i=1:Nsets
    NNBasis(:,:,i,i) = Vin(:,newInd,i);
    modelSpan(i) = numOfPoints;
    numOfPoints = numOfPoints + nPts;
end
modelSpan(Nsets+1) = numOfPoints;

CentroidPtsBelMod = zeros(numOfPoints, Nsets,'int32');


for i=1:Nsets
    
    spanI = modelSpan(i)+1:modelSpan(i+1);
    CentroidPtsBelMod(spanI,i) = 1:length(spanI);
    
    for j=i+1:Nsets
        
        spanJ = modelSpan(j)+1:modelSpan(j+1);
        % Pick only the correspondances that have the same neighboard
        % in both views
        
        Vi = Vin(:,newInd,i);
        Vj = Vin(:,newInd,j);
        
        % from i to j
        M = T(:, :, j, i);
        %R = closestRotation(M);
        Vc =  M * Vi;
        
        Corr1 = annquery(Vj, Vc, 1);
        Corr1 = newInd(Corr1);
        Corr1 = Corr1';
        %[Corr1, D1] = dsearchn(Model(j).vertices, Model(j).tri, Model(i).vertices);
        Corr1(:,2) = 1:length(Corr1);
        %            Corr1(D1 > tol,:) = [];
        
        % from j to i
        M = T(:, :, i, j);
        %R = closestRotation(M);
        Vc =  M * Vj;
        Corr2 = annquery(Vi, Vc, 1);
        Corr2 = newInd(Corr2);
        Corr2 = Corr2';
        %[Corr2, D2] = dsearchn(Model(i).vertices, Model(i).tri, Model(j).vertices);
        
        Corr2(:,2) = Corr2(:,1);
        Corr2(:,1) = 1:length(Corr2);
        %            Corr2(D2 > tol,:) = [];
        
        funique = ismember(Corr1,Corr2,'rows');
        Corr = Corr1(funique,:);
        
        CentroidPtsBelMod(spanI(Corr(:,2)),j) = Corr(:,1)';
        CentroidPtsBelMod(spanJ(Corr(:,1)),i) = Corr(:,2)';
        
        
    end
end

%CentroidPtsBelMod(sum(CentroidPtsBelMod>0,2)<2,:) = [];

for i=1:Nsets
    spanI = modelSpan(i)+1:modelSpan(i+1);
    tmp = CentroidPtsBelMod(spanI, :);
    tmp(sum(tmp>0,2)<Nsets,:) = [];
    for j=1:Nsets
        NNInd{j,i} = tmp(:,j);
    end
    
end



end