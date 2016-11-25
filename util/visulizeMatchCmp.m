% visualize sparse landmarks matching
function out = visulizeMatchCmp(dataList, consistentIdx1, consistentIdx2, BM1, BM2)

%BM1 = BM1(end:-1:1);
%BM2 = BM2(end:-1:1);
Nsets = length(dataList);
%Nsets = 2;
nLandmark = size(consistentIdx1, 1);
numEig = 30;
cyan_color = [1 1 1];
red_color = [1 0 1];

for i = 1:Nsets
    name = dataList{i};
    mesh = loadMeshLB( name, numEig);
    X{i} = mesh.vertices;
    T{i} = mesh.triangles;
    col{i} = repmat(cyan_color, mesh.nv, 1);
    %col{i} = X{i} - repmat(mean(X{i}), mesh.nv,1);
    %col{i} = col{i} ./ repmat(max(col{i}),mesh.nv,1);
end


maxX = max(X{1}(:,1));
minX = min(X{1}(:,1));
xshift = 2*(maxX - minX);

for i=2:Nsets
    X{i}(:,1) = X{i}(:,1) + (i-1)*xshift; 
    
end



h = figure;
% groudtruth
%h1 = subplot(1,2,1);
for i = 1:Nsets
    patch('Faces',T{i},'Vertices',X{i},'FaceColor','interp', ...
        'FaceVertexCData', col{i}, 'EdgeColor', 'none');
    axis equal; axis tight; axis off; cameratoolbar;

end

hold on
for i = 1:nLandmark
    pts1 = [];
    pts2 = [];
    pts3 = [];
    
    for j = 1:Nsets-1
        % PermSych
    si = consistentIdx1(BM1(i),j);
    ti = consistentIdx1(BM1(i),j+1);
    sp = X{j}(si, :);
    tp = X{j+1}(ti, :);
    pts1 = [pts1; sp; tp];
    % PermSych+
    si = consistentIdx2(BM2(i),j);
    ti = consistentIdx2(BM2(i),j+1);
    sp = X{j}(si, :);
    tp = X{j+1}(ti, :);
    pts2 = [pts2; sp; tp];
    % CSP
    %si = consistentIdx3(i,j);
    %ti = consistentIdx3(i,j+1);
    %sp = X{j}(si, :);
    %tp = X{j+1}(ti, :);
    %pts3 = [pts3; sp; tp];
    end
    
    %pts = [sp; tp];
    %ah=axes('position',[.2,.2,.6,.6],'visible','off');
    %h0 = text(0.5,0.9,['Best', num2str(i)], 'FontSize', 14);
    h1 = plot3(pts1(:,1), pts1(:,2), pts1(:, 3),'linewidth',2, 'LineStyle','--','Color','r','Marker','o');
    %l1 = legend('PermSych','Location','Best');
    %hold on
    %plot3(pts(:,1), pts(:,2), pts(:, 3),'linewidth',3, 'Color','w');
    %h2 = plot3(pts2(:,1), pts2(:,2), pts2(:, 3),'linewidth',2, 'LineStyle','-','Color','b','Marker','x');
    %l2 = legend('PermSych++','Location','Best');
    %plot3(pts(:,1), pts(:,2), pts(:, 3),'linewidth',3, 'Color','w');
    %h3 = plot3(pts3(:,1), pts3(:,2), pts3(:, 3),'linewidth',2, 'LineStyle','--','Color','g','Marker','s');
    %disp(consistentIdx1(BM1(i),:));

    %delete(h0);
    delete(h1);
    %delete(h2);
    %delete(h3);
    %h3 = text(pts2(:,1), pts2(:,2), pts2(:, 3), num2str(i));
end

% visualize landmarks
for i = 1:Nsets
    figure;
    
    col = repmat(cyan_color, mesh.nv, 1);
    for j=1:nLandmark
        col(consistentIdx2(j,i), :) = red_color;
    end
    
    patch('Faces',T{i},'Vertices',X{i},'FaceColor','interp', ...
        'FaceVertexCData', col, 'EdgeColor', 'none');
    axis equal; axis tight; axis off; cameratoolbar;

end

% for i=1:Nsets
%     P{i} = X{i}(consistentIdx2(:,i),:);
%     DT{i} = delaunayTriangulation(P{i});
%     faceColor  = [0.6875 0.8750 0.8984];
%     figure
%     tetramesh(DT{i},'FaceColor',faceColor,'FaceAlpha',0.3);
% end

out = 1;
end