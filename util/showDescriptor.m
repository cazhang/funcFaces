function h = showDescriptor(mesh, dloc, f, title)
if nargin < 4
    title = '';
end
X = mesh.vertices; T = mesh.triangles;
h = figure;
patch('vertices',X,'Faces',T,'FaceColor','interp','CData',f(:,dloc),'edgecolor','none'); 
axis equal; cameratoolbar; hold on; colorbar
set(h,'WindowStyle','docked','name',title);    
colormap;
