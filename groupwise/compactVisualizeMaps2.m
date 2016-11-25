% visualize maps for a pair of shapes
function compactVisualizeMaps2(mesh1, basis1, mesh2, basis2, Fmaps)

h1s = [];
h2s = [];
h3s = [];
h4s = [];
h5s = [];
h6s = [];


for i=1:length(Fmaps)
    [h1,h2,h3,h4,h5,h6] = compactVisualizeMap2(mesh1, basis1, mesh2, basis2, Fmaps{i});
    h1s = [h1s, h1]; 
    h2s = [h2s, h2];
    h3s = [h3s, h3];
    h4s = [h4s, h4];
    h5s = [h5s, h5];
    h6s = [h6s, h6];
end

for i=1:length(Fmaps)
    hlink = linkprop(h1s(i),{'CameraPosition','CameraUpVector'});
    key = 'graphics_linkprop';
    setappdata(h1s(1),key,hlink); 
    
    hlink = linkprop(h2s(i),{'CameraPosition','CameraUpVector'});
    key = 'graphics_linkprop';
    setappdata(h2s(1),key,hlink);  
    
    hlink = linkprop(h3s(i),{'CameraPosition','CameraUpVector'});
    key = 'graphics_linkprop';
    setappdata(h3s(1),key,hlink);  
    
    hlink = linkprop(h4s(i),{'CameraPosition','CameraUpVector'});
    key = 'graphics_linkprop';
    setappdata(h4s(1),key,hlink);  
    
    hlink = linkprop(h5s(i),{'CameraPosition','CameraUpVector'});
    key = 'graphics_linkprop';
    setappdata(h5s(1),key,hlink);  
    
    hlink = linkprop(h6s(i),{'CameraPosition','CameraUpVector'});
    key = 'graphics_linkprop';
    setappdata(h6s(1),key,hlink);  
    
    
end