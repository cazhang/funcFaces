function matrix = qrMap(F, G, epsilon)
% Computes a map matrix using QR decomposition.

[q1 r1] = qr(F);
[q2 r2] = qr(G);

r3 = r2*pinv(r1,epsilon);
[u1 s1 v1] = svd(r3);
d1 = diag(s1);
d1(d1>epsilon) = 1;

matrix =  q2*u1*diag(d1)*v1'*q1';