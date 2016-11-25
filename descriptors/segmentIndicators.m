function F = segmentIndicators(mesh, segs);
F = [];
if isempty(segs)
    return;
end
for i=1:length(segs)
    C = (mesh.segmentation == segs(i));
    if mesh.segMethod.faceSeg
        C = faceToVertex(C,mesh);
    end
    F = [F C];
end

