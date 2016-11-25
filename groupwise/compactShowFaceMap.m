function [h, h1,h2,h3,h4,h5,h6] = compactShowFaceMap(mesh1, basis1, mesh2, basis2, fmap, dest)
h = figure; 
X1 = mesh1.vertices; T1 = mesh1.triangles;
X2 = mesh2.vertices; T2 = mesh2.triangles;

tex1 = mesh1.texture;  
tex1 = tex1./255;

tex2 = mesh2.texture;
tex2 = tex2./255;

% groundtruth tex1
h1 = subplot(1,6,1);
patch('Faces',T1,'Vertices',X1,'FaceColor','interp', ...
      'FaceVertexCData', tex1, 'EdgeColor', 'none'); 
axis equal; axis tight; axis off; cameratoolbar; 

% groundtruth tex2
h2 = subplot(1,6,2);
%col2 = transferFunction(col1, basis1, basis2, fmap);
patch('Faces',T2,'Vertices',X2,'FaceColor','interp', ...
      'FaceVertexCData', tex2, 'EdgeColor', 'none'); 
axis equal; axis tight; axis off; cameratoolbar;          

% pairwise funcmap
h3 = subplot(1,6,3);
%col2 = transferFunction(col1, basis1, basis2, fmap);
patch('Faces',T1,'Vertices',X1,'FaceColor','interp', ...
      'FaceVertexCData', tex2(dest{1},:), 'EdgeColor', 'none'); 
axis equal; axis tight; axis off; cameratoolbar;         

% groupwise funcmap
h4 = subplot(1,6,4);
%col2 = transferFunction(col1, basis1, basis2, fmap);
%col4 = col1(dest{3},:);
patch('Faces',T1,'Vertices',X1,'FaceColor','interp', ...
      'FaceVertexCData', tex2(dest{2},:), 'EdgeColor', 'none'); 
axis equal; axis tight; axis off; cameratoolbar;          

% pairwise icp
h5 = subplot(1,6,5);
%col2 = transferFunction(col1, basis1, basis2, fmap);
%col5 = col1(dest{4},:);
patch('Faces',T1,'Vertices',X1,'FaceColor','interp', ...
      'FaceVertexCData', tex2(dest{3},:), 'EdgeColor', 'none'); 
axis equal; axis tight; axis off; cameratoolbar;  

% groupwise icp
h6 = subplot(1,6,6);
%col2 = transferFunction(col1, basis1, basis2, fmap);
%col5 = col1(dest{4},:);
patch('Faces',T1,'Vertices',X1,'FaceColor','interp', ...
      'FaceVertexCData', tex2(dest{4},:), 'EdgeColor', 'none'); 
axis equal; axis tight; axis off; cameratoolbar;  
  



set(h,'WindowStyle','normal');  




