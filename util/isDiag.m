function d = isDiag(M)
if norm(M - diag(diag(M))) < 1e-10
    d = 1;
else
    d = 0;
end
    
