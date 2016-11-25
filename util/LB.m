function L = LB(mesh)
A = spdiags(1./mesh.areaWeights, 0, mesh.nv, mesh.nv);
L = A*mesh.cotLaplace;