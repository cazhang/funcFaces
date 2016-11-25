function Cout = icpFmap(Cin, basis1, basis2, step, numIter)

if nargin < 5
    numIter = 10;
end
if nargin < 4
    step = 1;
end
   

Cout = closestRotation(Cin);

V1 = basis1;
V2 = basis2;

testidx = 1:step:size(V1,1);
V1 = V1(testidx,:);
V2 = V2(:,:);

for k = 1:numIter
    Vc = Cout*V1';
    fprintf('ann ... ');
    t = tic;
    nnidx = annquery(V2',Vc, 1);
    fprintf('%g seconds\n',toc(t));

    W = V2(nnidx,:)'*V1(:,:);
    [uu ss vv] = svd(W);
    Cout = uu*vv';
end
