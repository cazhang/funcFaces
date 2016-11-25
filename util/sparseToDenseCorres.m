% Based on sparse landmark membership, further assign denser level
function [DenseConsistentVtx] = sparseToDenseCorres(Nlandmark,V,ConsistentVtx,SAMPLES,MEMBERS)
lm_set = [];
Num_set = [];% shape x landmarks


for area_id = 1:Nlandmark
    area_SAMPLES = [];
    lm_num = ConsistentVtx( area_id, :);
    for shape_id = 1:Nsets
        lm_id = find(SAMPLES(:,shape_id) == lm_num(shape_id));
        memb = MEMBERS(:, shape_id);
        lm_set{shape_id, area_id} = find(memb == lm_id);
        Num_set(shape_id, area_id) = length(lm_set{shape_id, area_id});
    end
    
    % choose subset of basis
    num_sel = min(Num_set(:, area_id));
    for shape_id = 1:Nsets
        tmp = randperm(Num_set(shape_id, area_id));
        samples = lm_set{shape_id, area_id}(tmp(1:num_sel));
        area_SAMPLES = [area_SAMPLES samples];
        Vsel(:,:,shape_id) = V(:,samples,shape_id);
        clear tmp
    end
    % run PermSych
    Nlandmark = num_sel;
    [evnoise, recPM, gs] = permSych(Nsets,Nlandmark,T,Vsel);
    AreaConsistentVtx = area_SAMPLES(:,1);
    for i=1
        for j=i+1:Nsets
            tmp = area_SAMPLES(:,j);
            tmp = tmp(gs{i}{j});
            AreaConsistentVtx = [AreaConsistentVtx tmp];
        end
    end
    DenseConsistentVtx{area_id} = AreaConsistentVtx;
    %out = visulizeMatch(dataList, AreaConsistentVtx);
    
    
end
