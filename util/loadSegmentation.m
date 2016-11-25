function mesh = loadSegmentation(mesh, segMethod)
filename = [mesh.name '_' segMethod.name '_seg.txt'];
seg = dlmread(filename);
%load(filename);
mesh.segMethod = segMethod;
if min(seg) == 0
    seg = seg+1;
    dlmwrite(filename,seg,'newline','pc');
end
mesh.segmentation = seg;
    