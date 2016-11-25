function EL = edgeLengthMatrix(mesh)
nv = mesh.nv;
X = mesh.vertices;
T = mesh.triangles;

I = [T(:,1);T(:,2);T(:,3)];
J = [T(:,2);T(:,3);T(:,1)];
S = normv(X(I,:) - X(J,:));
EL = sparse(double(I), double(J), double(S), nv,nv);

function nn = normv(V)
nn = sqrt(sum(V.^2,2));