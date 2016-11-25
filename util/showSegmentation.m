function h = showSegmentation(mesh, segMethod)
X = mesh.vertices; T = mesh.triangles;
segmentation = dlmread([mesh.name '_' segMethod.name '_seg.txt']);
h = figure;
if segMethod.faceSeg 
    faceColor = 'flat';
else
    faceColor = 'interp';
end
patch('vertices',X,'Faces',T,'FaceColor',faceColor,'CData',segmentation,'edgecolor','none'); 
axis equal; cameratoolbar; hold on; colorbar
set(h,'WindowStyle','docked');    
colormap(segColormap(max(segmentation)));

ns = max(segmentation);
ticks = 1+[0:ns-1]*ns/(ns+1);
diff = ticks(2)-ticks(1);
ticks = ticks + diff/2;
tickLabels = cellfun(@num2str, num2cell([1:ns]),'UniformOutput',false);

h = colorbar;
set(h,'Ytick',ticks,'YtickLabel',tickLabels);


