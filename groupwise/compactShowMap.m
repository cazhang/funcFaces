function [h, h1,h2,h3,h4,h5,h6] = compactShowMap(mesh1, basis1, mesh2, basis2, fmap, dest)
h = figure; 
X1 = mesh1.vertices; T1 = mesh1.triangles;
X2 = mesh2.vertices; T2 = mesh2.triangles;
col1 = X1 - repmat(mean(X1),mesh1.nv,1);    
col1 = col1./repmat(max(col1),mesh1.nv,1);

% groudtruth
h1 = subplot(1,6,1);
patch('Faces',T1,'Vertices',X1,'FaceColor','interp', ...
      'FaceVertexCData', col1, 'EdgeColor', 'none'); 
axis equal; axis tight; axis off; cameratoolbar; 

% pairwise funcmap
h2 = subplot(1,6,2);
%col2 = transferFunction(col1, basis1, basis2, fmap);
col2 = col1(dest{1},:);
patch('Faces',T2,'Vertices',X2,'FaceColor','interp', ...
      'FaceVertexCData', col2, 'EdgeColor', 'none'); 
axis equal; axis tight; axis off; cameratoolbar;          

% groupwise funcmap
h3 = subplot(1,6,3);
%col2 = transferFunction(col1, basis1, basis2, fmap);
col3 = col1(dest{2},:);
patch('Faces',T2,'Vertices',X2,'FaceColor','interp', ...
      'FaceVertexCData', col3, 'EdgeColor', 'none'); 
axis equal; axis tight; axis off; cameratoolbar;         

% groupwise icp
h4 = subplot(1,6,4);
%col2 = transferFunction(col1, basis1, basis2, fmap);
col4 = col1(dest{3},:);
patch('Faces',T2,'Vertices',X2,'FaceColor','interp', ...
      'FaceVertexCData', col4, 'EdgeColor', 'none'); 
axis equal; axis tight; axis off; cameratoolbar;          

% pairwise icp
h5 = subplot(1,6,5);
%col2 = transferFunction(col1, basis1, basis2, fmap);
col5 = col1(dest{4},:);
patch('Faces',T2,'Vertices',X2,'FaceColor','interp', ...
      'FaceVertexCData', col5, 'EdgeColor', 'none'); 
axis equal; axis tight; axis off; cameratoolbar;    

%h6 = subplot(1,6,6);
%col2 = transferFunction(col1, basis1, basis2, fmap);
%col6 = col1(dest{5},:);
%patch('Faces',T2,'Vertices',X2,'FaceColor','interp', ...
%      'FaceVertexCData', col6, 'EdgeColor', 'none'); 
%axis equal; axis tight; axis off; cameratoolbar;    



set(h,'WindowStyle','normal');  




