function S = symmetryOperatorGT(mesh, basis)
s=load([mesh.name '_sym']);
Ss = sparse(1:mesh.nv,s.corr,ones(mesh.nv,1),mesh.nv,mesh.nv);
S = descriptorsToCoeffs(Ss*basis, basis);
S(abs(S)<1e-3)=0;
