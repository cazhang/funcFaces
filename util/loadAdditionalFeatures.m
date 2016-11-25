function mesh = loadAdditionalFeatures(mesh)

fprintf(1,'Loading additional features...\n');
mesh.additionalFeatures = dlmread([mesh.name '_features.txt'],' ');