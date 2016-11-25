function matrix = l1RegularizedMap(F, G, lambda)
% Computes a map matrix by solving:
%   min .5*||AF-G||_2^2 + lambda * ||A||_1

A = F';
B = G';
numCols = size(B,2);

variableMatrix = sparse(1:numCols,1:numCols,repmat(1,1,numCols),numCols,2*numCols);
absMatrix = sparse(1:numCols,(1:numCols)+numCols,repmat(1,1,numCols),numCols,2*numCols);

H = variableMatrix' * A' * A * variableMatrix;

matrixIneq = [-absMatrix+variableMatrix ; -absMatrix-variableMatrix];
rhsIneq = zeros(2*numCols,1);

X = zeros(size(A,2), size(B,2));

for i=1:numCols
    bCol = B(:,i);
    
    f = zeros(2*numCols,1);
    f(1:numCols,1) = -A'*bCol;
    f((numCols+1):(2*numCols)) = lambda;
    
    Xcol = cplexqcp(H,f,matrixIneq,rhsIneq);
    column = Xcol(1:numCols);
    
    X(:,i) = column;
end
matrix = X';