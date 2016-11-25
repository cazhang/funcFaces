function [err, bestMatch] = quantifyMatch(dataList, consistentIdx, V, T)
Nsets = length(dataList);
nLandmark = size(consistentIdx, 1);
numEig = 30;
err = 0;
cyan_color = [0 1 1];
accDistance = zeros(nLandmark, 1);
for i = 1:Nsets
    for j=1:Nsets
        if i ~= j
%             name2 = dataList{j};
%             mesh2 = loadMeshLB( name2, numEig);
%             X2 = mesh2.vertices;
%             edgeLengths2 = mesh2.edgeLens;
%             meanEdgeLength2 = mean(nonzeros(edgeLengths2));
%             
%             toTest = consistentIdx(:, i);
%             destinations = consistentIdx(:, j);
%             distances = normv(X2(toTest,:) - X2(destinations,:));
                

            V1 = V(:,consistentIdx(:,i),i);
            V2 = V(:,consistentIdx(:,j),j);
            Vproj = T(:,:,j,i)*V1;
            distances = normv(V2' - Vproj');

            %normDistance = distances / meanEdgeLength2;
            accDistance = accDistance + distances;
            err = err + mean(distances);
        end
    end
    %col{i} = X{i} - repmat(mean(X{i}), mesh.nv,1);
    %col{i} = col{i} ./ repmat(max(col{i}),mesh.nv,1);
end

[junk, bestMatch] = sort(accDistance);