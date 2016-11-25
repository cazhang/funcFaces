function vertexSignal = faceToVertex(faceSignal, mesh)
% Converts a signal on faces to a signal on vertices on a mesh.
% Averages the signal of near by faces, weighted by the face area
% 

f = faceSignal;
X = mesh.vertices;
T = mesh.triangles;

normals = cross(X(T(:,1),:) - X(T(:,2),:), X(T(:,1),:) - X(T(:,3),:));
areas = sqrt(sum(normals.^2,2))/2;

nf = size(T,1);
nv = size(X,1);

I = [T(:,1);T(:,2);T(:,3)];
J = [1:nf,1:nf,1:nf]';
S = areas(J);
f2v_w = sparse(I,J,S,nv,nf);
f2v_w = spdiags(1./sum(f2v_w,2),0,nv,nv)*f2v_w;

vertexSignal = f2v_w*f;