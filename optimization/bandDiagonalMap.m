function matrix = bandDiagonalMap(F, G, band)
% Computes a map matrix restricted to be in a diagonal of the given width.

N = size(F,1);
matrix = zeros(N);

% True bandsize = band*2+1
for i=1:N
    start = i - band;
    stop = i + band;

    start = max(start, 1);
    stop = min(stop, N);
    
    lazy = F(start:stop,:)'\G(i,:)';
    matrix(i,:) = [zeros(1,start-1) lazy' zeros(1, N-stop)];
end
