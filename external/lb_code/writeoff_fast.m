function [T X Y Z] = writeoff_fast(filename,T,X,Y,Z)
    fid = fopen(filename,'wt');

    if fid==-1 
        error('File not found or permission denied');
    end

    if size(X,2) ~= 1 || size(X,1) ~= size(Y,1)
        error('dimensions mismatch');
    end

    fprintf(fid,'OFF\n');

fprintf(fid,'%d %d 0\n',size(X,1),size(T,1));
fclose(fid);

C = [X Y Z];
fprintf('writing coordinates\n');
save(filename, 'C', '-append', '-ascii','-double');
%dlmwrite(filename,C,'-append','delimiter','\t','precision','%.9f');

fprintf('writing triangles\n');
C = [3*ones(size(T,1),1) T-1];
dlmwrite(filename,C,'-append','delimiter','\t');

end
