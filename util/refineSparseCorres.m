% use PermSych to recover one-to-one point correspondence
function [Corres, ConsistentVtx, refinedConsistentVtx, err,SAMPLES,MEMBERS, bestMatchPerm, bestMatchRefine] = refineSparseCorres(classID, dataList,X, Nlandmark, V, Nsets)

T = M2T(X);
flag = 0;

samples_path = '/local/sda7/cz679/funMap/CSP_Codes/samples/';
SAMPLES = [];
MEMBERS = [];
Vsel = [];
for i = 1:Nsets
    load( [samples_path, 'FAUST', num2str(classID), '_', num2str(i),'_rnd',num2str(Nlandmark),'.mat']);
    SAMPLES = [SAMPLES samples'];
    MEMBERS = [MEMBERS members'];
    Vsel(:,:,i) = V(:, samples, i);
end

% fine the correspondence in the remaining shapes
Corres = SAMPLES(:,1);
for i=1:Nsets-1
    toTest = Corres(:,i)';
    basis1 = V(:,:,i)';
    basis2 = V(:,:,i+1)';
    fmap = T(:,:,i+1,i);
    dest = mapClosestDelta(fmap, toTest, basis1, basis2);
    Corres = [Corres dest'];
end

% find the group consistent assignment of landmarks
[evnoise, recPM, gs] = permSych(Nsets,Nlandmark,T,Vsel);

ConsistentVtx = SAMPLES(:,1);
for i=1
    for j=i+1:Nsets
        tmp = SAMPLES(:,j);
        tmp = tmp(gs{i}{j});
        ConsistentVtx = [ConsistentVtx tmp];
    end
end
clear Vsel

% % remove false correspondence
% reducedConsis = [];
% for i=1:Nlandmark
%     membership = MEMBERS(ConsistentVtx(i,1),1);
%     corres = ConsistentVtx(i,:);
%     corres_member = MEMBERS( corres, 1);
%     if unique(corres_member) == membership
%         reducedConsis = [reducedConsis; ConsistentVtx(i,:)];
%     end
% end
% refinedConsistentVtx = reducedConsis;

% visualize
%out = visulizeMatch(dataList, ConsistentVtx);
% quantify error
err = [];
[err(1), bestMatchPerm] = quantifyMatch(dataList, ConsistentVtx, V, T);
bestMatchRefine = bestMatchPerm;

refinedConsistentVtx = ConsistentVtx(1:Nlandmark/2,:);

while flag
%     occp = [];
%     for k=2:Nsets
%         occp{k} = zeros(Nlandmark,1);
%     end
    % adjust assignment by comparing with mean value in embeded space
    for i=1:Nlandmark
        meanV = 0;
        for j=1:Nsets
            % need to project with map to shape 1
            meanV = meanV + T(:,:,1,j)*V(:,refinedConsistentVtx(i,j),j);
        end
        meanV = meanV ./ Nsets;
        %meanV = V(:,refinedConsistentVtx(i, 1),1);
        %newConsistentVtx(:,1) = refinedConsistentVtx(:,1);
        for j=1:Nsets
            if flag==1
                % find in the landmark set
                dist = pdist2(V(:,refinedConsistentVtx(:,j),j)'*T(:,:,1,j)', meanV');
            elseif flag==2
                % find in the whole set
                dist = pdist2(V(:,:,j)'*T(:,:,1,j)', meanV');
            end
            %dist(find(occp{j}>0)) = inf;
            [a,b] = min(dist);
            %occp{j}(b,1) = 1;
            if b == i
                newConsistentVtx(i,j) = refinedConsistentVtx(i,j);
            else
                if flag==1
                    newConsistentVtx(i,j) = refinedConsistentVtx(b,j);
                    %newConsistentVtx(b,j) = ConsistentVtx(i,j);
                elseif flag==2
                    newConsistentVtx(i,j) = b;
                end
            end
        end
    end
    refinedConsistentVtx = newConsistentVtx;
    
    [err(end+1), bestMatchRefine] = quantifyMatch(dataList, refinedConsistentVtx, V, T);
    
    if err(end-1)-err(end) < 1e-6
        break;
    end
end

% post-processing (remove redundent correspondence)

