%% evaluate mapping, generate CMC curve

clc;
clear all;
close all;

setupFunctionalMaps
Nsets = 11;
numEig = 30;
dataList = cell(Nsets);
for i=1:Nsets
    dataList{i} = ['cat',num2str(i-1)];
end

Dist1 = [];
Dist2 = [];
Dist3 = [];
Dist4 = [];

for i = 1:Nsets
    for j = 1:Nsets
        if i~=j && i < j
            name2 = dataList{i};
            name1 = dataList{j};
            p2p_file_name = [name1, '_to_', name2, '.mat'];
            load( p2p_file_name );
            
            mesh1 = loadMeshLB( name1, numEig);
            mesh2 = loadMeshLB( name2, numEig);
            
            edgeLengths2 = mesh2.edgeLens;
            meanEdgeLength2 = mean(nonzeros(edgeLengths2));
            X2 = mesh2.vertices;
            distances1 = normv(X2(:,:) - X2(destinations{1},:));
            distances2 = normv(X2(:,:) - X2(destinations{2},:));
            distances3 = normv(X2(:,:) - X2(destinations{3},:));
            distances4 = normv(X2(:,:) - X2(destinations{4},:));
            
            distances1 = distances1 / meanEdgeLength2;
            distances2 = distances2 / meanEdgeLength2;
            distances3 = distances3 / meanEdgeLength2;
            distances4 = distances4 / meanEdgeLength2;
            
            Dist1 = [Dist1 distances1];
            Dist2 = [Dist2 distances2];
            Dist3 = [Dist3 distances3];
            Dist4 = [Dist4 distances4];
            
            
        end
    end
end

meanDis1 = mean(Dist1, 2);
meanDis2 = mean(Dist2, 2);
meanDis3 = mean(Dist3, 2);
meanDis4 = mean(Dist4, 2);

percent1 = [];
percent2 = [];
percent3 = [];
percent4 = [];
for i=1:100
    a = find(meanDis1 < i*0.1);
    percent1(i) = length(a) / length(meanDis1);
    a = find(meanDis2 < i*0.1);
    percent2(i) = length(a) / length(meanDis2);
    a = find(meanDis3 < i*0.1);
    percent3(i) = length(a) / length(meanDis3);
    a = find(meanDis4 < i*0.1);
    percent4(i) = length(a) / length(meanDis4);
    
end


figure;

plot(percent1, 'r');
hold on
plot(percent2, 'b');
hold on
plot(percent3, 'r.');
hold on
plot(percent4, 'b.');

      