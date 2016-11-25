% use PermSych to recover one-to-one point correspondence
function [ConsistentVtx, refinedConsistentVtx, err] = refineCSPCorres(classID, dataList,X, Nlandmark, V, Nsets)

T = M2T(X);
flag = 2;

% load CSP
CSP_data = load(['FAUST_Corre_CSP_n', num2str(Nlandmark),'_',num2str(classID),'.mat']);
ConsistentVtx = CSP_data.consistentVertexIds;
% visualize
%out = visulizeMatch(dataList, ConsistentVtx);
% quantify error
err = [];
err(1) = quantifyMatch(dataList, ConsistentVtx);

refinedConsistentVtx = ConsistentVtx;
while 1
    % adjust assignment by comparing with mean value in embeded space
    for i=1:Nlandmark
        refRnd = randperm(Nsets);
        refId = 1;
        
        meanV = 0;
        meanV = V(:,refinedConsistentVtx(i, 1),1);
        %for j=1:Nsets
            % need to project with map to shape 1
        %    meanV = meanV + T(:,:,refId,j)*V(:,refinedConsistentVtx(i,j),j);
        %end
        %meanV = meanV ./ Nsets;
        newConsistentVtx(:,1) = refinedConsistentVtx(:,1);
        for j=2:Nsets
            if flag==1
                % find in the landmark set
                dist = pdist2(V(:,refinedConsistentVtx(:,j),j)'*T(:,:,refId,j)', meanV');
            elseif flag==2
                % find in the whole set
                dist = pdist2(V(:,:,j)'*T(:,:,refId,j)', meanV');
            end
            [a,b] = min(dist);
            
            if flag==1
                newConsistentVtx(i,j) = refinedConsistentVtx(b,j);
                newConsistentVtx(b,j) = refinedConsistentVtx(i,j);
            elseif flag==2
                newConsistentVtx(i,j) = b;
            end
            
        end
    end
    refinedConsistentVtx = newConsistentVtx;
    
    err(end+1) = quantifyMatch(dataList, refinedConsistentVtx);
    % terminate loop
    if err(end-1)-err(end) < 1e-6
        break;
    end
end