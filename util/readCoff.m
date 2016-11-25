function [X,T,H] = readCoff(filename, m)

fid = fopen(filename,'r');
if( fid==-1 )
    error('Cannot open the file.');
    return;
end

str = fgets(fid);   % -1 if eof
if ~strcmp(str(1:4), 'COFF')
    error('The file is not a valid COFF one.');    
end

str = fgets(fid);
sizes = sscanf(str, '%d %d %d\n', 3);
nv = sizes(1);
nf = sizes(2);

% Read vertices
str = '%lf %lf %lf ';
for i=1:m
    str = [str '%lf '];
end
str = [str '\n'];
[X,cnt] = fscanf(fid, str, [3+m,nv]);
if cnt~=(3+m)*nv
    warning('Problem in reading vertices.');
end
X = X';

H = X(:,4:end);
X = X(:,1:3);

[T,cnt] = fscanf(fid,'3 %ld %ld %ld\n', [3,inf]);
T = T'+1;

fclose(fid);