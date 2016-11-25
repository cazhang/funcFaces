function X = initEye(d, N)

X = zeros(d, d, N-1);

for i=1:N-1
    X(:,:,i) = eye(d);
end




end