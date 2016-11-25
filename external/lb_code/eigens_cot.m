% Use this file to compute the _cot.mat files for given shape files.

nev = 300;

a = dir('../shapes/*.mat');

addpath('./MeshLP');

for i=1:length(a)
    prefix = a(i).name;
	prefix = prefix(1:end-4);
       
    fprintf('%s\n',prefix)
    outfile = sprintf('../shapes/%s_cot.mat', prefix);
    if(~exist(outfile,'file') && ~endswith(prefix,'_cot'))
        fprintf('%s does not exist, computing',outfile);
        filename=sprintf('../shapes/%s.mat', prefix);
        
        A = load(filename);
        
        offfilename=sprintf('../shapes/%s.off', prefix);
        
        A = A.surface;
        
        writeoff_fast(offfilename, A.TRIV, A.X, A.Y, A.Z);
        
        filename = offfilename;
        
        [W A]=cotlp_matrix(filename);
        
        Am = sparse([1:length(A)], [1:length(A)], A);
        
        [evecs evals] = eigs(W, Am, nev, -1e-5);
        evals = diag(evals);
        
        filename=sprintf('./%s_cot', prefix);
        save(filename, 'W', 'A', 'evecs', 'evals');
      end
end
