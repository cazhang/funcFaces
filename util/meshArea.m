function area = meshArea(mesh)

X = mesh.vertices; T = mesh.triangles;
Ar = normv(cross(X(T(:,1),:) - X(T(:,2),:), X(T(:,1),:) - X(T(:,3),:)))/2;
area = sum(Ar);