function [out] = displayFaceResult(Fmaps, x0, xsol, xicp, Nsets, Ndims, nICP, isShow, toSave, dataList, nn_method)
% pairwise
T0 = M2T( x0 );
% groupwise
That = cell(length(xsol),1);
for i = 1:length(xsol)
    That{i} = M2T( xsol{i});
end
%Ticp = M2T( xicp{end} );
% group ICP
Ticp = cell(length(xicp),1);
for i=1:length(xicp)
    Ticp{i} = M2T( xicp{i} );
end


%% display before
%fprintf('Objective value: %.6e\n', O);

% figure;
% for i = 1: Nsets
%     for j = 1: Nsets
%         subplot(Nsets, Nsets, i + (j - 1) * Nsets)
%         %tmp = (T(:, :, i, j) - That(:, :, i, j)).^2;
%         tmp = T0(:, :, i, j);
%         imagesc(tmp), colorbar
%         tmp = sum(tmp(:)) / (Ndims * Ndims);
%         if tmp == 0
%             caxis([0 1]);
%         end
%         title('After pairwise');
%         %title(sprintf('That_{%d%d} vs T_{%d%d} (%.2e)', i, j, i, j, tmp));
%         %fprintf('Estimation error %d-%d: %.2e\n', i, j, tmp);
%         set(gca, 'XTick', [])
%         set(gca, 'YTick', [])
%     end
% end

%         for i = 1: Nsets
%             for j = 1: Nsets
%                 if i ~= j
%                     tmp = (T0(:, :, i, j) * T0(:, :, j, i) - eye(Ndims)).^2;
%                     fprintf('Invertibility constraint %d-%d: %.2e\n', i, j, sum(tmp(:)) / (Ndims * Ndims))
%                 end
%             end
%         end
%
%         for i = 1: Nsets
%             for j = 1: Nsets
%                 for k = 1: Nsets
%                     if i ~= j && j ~= k && i ~= k
%                         tmp = (T0(:, :, i, j) * T0(:, :, j, k) - T0(:, :, i, k)).^2;
%                         fprintf('Transitivity constraint %xicpd-%d-%d: %.2e\n', i, j, k, sum(tmp(:)) / (Ndims * Ndims))
%                     end
%                 end
%             end
%         end


%% display after


% figure;
% for i = 1: Nsets
%     for j = 1: Nsets
%         subplot(Nsets, Nsets, i + (j - 1) * Nsets)
%         %tmp = (T(:, :, i, j) - That(:, :, i, j)).^2;
%         tmp = That(:, :, i, j);
%         imagesc(tmp), colorbar
%         tmp = sum(tmp(:)) / (Ndims * Ndims);
%         if tmp == 0
%             caxis([0 1]);
%         end
%         title('After groupwise');
%         %title(sprintf('That_{%d%d} vs T_{%d%d} (%.2e)', i, j, i, j, tmp));
%         %fprintf('EstimationO error %d-%d: %.2e\n', i, j, tmp);
%         set(gca, 'XTick', [])
%         set(gca, 'YTick', [])
%     end
% end

% for i = 1: Nsets
%     for j = 1: Nsets
%         if i ~= j
%             tmp = (That(:, :, i, j) * That(:, :, j, i) - eye(Ndims)).^2;
%             fprintf('Invertibility constraint %d-%d: %.2e\n', i, j, sum(tmp(:)) / (Ndims * Ndims))
%         end
%     end
% end
% 
% for i = 1: Nsets
%     for j = 1: Nsets
%         for k = 1: Nsets
%             if i ~= j && j ~= k && i ~= k
%                 tmp = (That(:, :, i, j) * That(:, :, j, k) - That(:, :, i, k)).^2;
%                 fprintf('Transitivity constraint %d-%d-%d: %.2e\n', i, j, k, sum(tmp(:)) / (Ndims * Ndims))
%             end
%         end
%     end
% end
% 
% clear i j k tmp O

% for k=1:length(xicp)
% figure;
% for i = 1: Nsets
%     for j = 1: Nsets
%         subplot(Nsets, Nsets, i + (j - 1) * Nsets)
%         %tmp = (T(:, :, i, j) - That(:, :, i, j)).^2;
%         tmp = Ticp{k}(:, :, i, j);
%         imagesc(tmp), colorbar
%         tmp = sum(tmp(:)) / (Ndims * Ndims);
%         if tmp == 0
%             caxis([0 1]);
%         end
%         title(sprintf('After ICP, %d', k));
%         %title(sprintf('That_{%d%d} vs T_{%d%d} (%.2e)', i, j, i, j, tmp));
%         %fprintf('EstimationO error %d-%d: %.2e\n', i, j, tmp);
%         set(gca, 'XTick', [])
%         set(gca, 'YTick', [])
%     end
% end
% end

%% compare pairwise maps with groupwise maps
% only compare i->1 mapping

% % set mesh and basis
% for i=2:Nsets
%     for j=1:Nsets
%         if i~=j % dont bother diagnal elements
%
%             if j==1
%                 Fmaps{i,1}{1,1}.mesh1 = Fmaps{1,i}{1,1}.mesh2;
%                 Fmaps{i,1}{1,1}.mesh2 = Fmaps{1,i}{1,1}.mesh1;
%                 Fmaps{i,1}{1,1}.basis1 = Fmaps{1,i}{1,1}.basis2;
%                 Fmaps{i,1}{1,1}.basis2 = Fmaps{1,i}{1,1}.basis1;
%             else
%                 Fmaps{i,j}{1,1}.mesh1 = Fmaps{1,j}{1,1}.mesh1;
%                 Fmaps{i,j}{1,1}.mesh2 = Fmaps{1,i}{1,1}.mesh1;
%                 Fmaps{i,j}{1,1}.basis1 = Fmaps{1,j}{1,1}.basis1;
%                 Fmaps{i,j}{1,1}.basis2 = Fmaps{1,i}{1,1}.basis1;
%
%             end
%
%         end
%     end
% end

for i=1:Nsets
    for j=1:Nsets
        if i ~= j
        % set initial map (only first row pairwise init)
        %Fmaps{i,j}{1,1}.maps{1} = T0(:,:,i,j);
        Fmaps{i,j}{1,1}.params{1} = 0;
        %Fmaps{i,j}{1,1}.method{1} = 'PairFunc';
        % set groupwise map
        Fmaps{i,j}{1,1}.maps{2} = That{1}(:,:,i,j);
        Fmaps{i,j}{1,1}.params{2} = 1;
        %Fmaps{i,j}{1,1}.method{2} = 'GroupFunc';
        
        
        %Fmaps{i,j}{1,1}.maps{3} = That{2}(:,:,i,j);
        %Fmaps{i,j}{1,1}.params{3} = 1;
%         Fmaps{i,j}{1,1}.method{3} = 'GroupId';
%         
%         Fmaps{i,j}{1,1}.maps{4} = That{3}(:,:,i,j);
%         Fmaps{i,j}{1,1}.params{4} = 1;
%         Fmaps{i,j}{1,1}.method{4} = 'GroupBIM';
        
        
        % set groupicp map
        %Fmaps{i,j}{1,1}.maps{end+1} = Ticp(:,:,i,j);
        %Fmaps{i,j}{1,1}.params{end+1} = 2;
        %Fmaps{i,j}{1,1}.method{end+1} = 'GroupICP';
        %         for k = 1:length(Ticp)
        %             Fmaps{i,j}{1,1}.maps{end+1} = Ticp{k}(:,:,i,j);
        %             Fmaps{i,j}{1,1}.params{end+1} = 2;
        %         end
        
        
        % set pairwise icp map
        if isfield(Fmaps{i,j}{1,1},'ppMaps')
            Fmaps{i,j}{1,1}.maps{3} = Fmaps{i,j}{1,1}.ppMaps{1};
            Fmaps{i,j}{1,1}.params{3} = 2;
            %Fmaps{i,j}{1,1}.method{3}= 'PairICP';
            Fmaps{i,j}{1,1} = rmfield(Fmaps{i,j}{1,1}, 'ppMaps');
        else
            Fmaps{i,j}{1,1}.maps{3} = Fmaps{i,j}{1,1}.maps{1};
            Fmaps{i,j}{1,1}.params{3} = 2;
            
        end
        
        Fmaps{i,j}{1,1}.maps{4} = Ticp{1}(:,:,i,j);
        Fmaps{i,j}{1,1}.params{4} = 3;
        %Fmaps{i,j}{1,1}.maps{5} = Ticp{2}(:,:,i,j);
        %Fmaps{i,j}{1,1}.params{5} = 3;
        %Fmaps{i,j}{1,1}.maps{6} = Ticp{3}(:,:,i,j);
        %Fmaps{i,j}{1,1}.params{6} = 3;
        %Fmaps{i,j}{1,1}.maps{7} = Ticp{4}(:,:,i,j);
        %Fmaps{i,j}{1,1}.params{7} = 3;
        
        %Fmaps{i,j}{1,1}.method{4}= 'GroupICP';
        end
    end
end

plist_fname = '/benchmark/PairList/VerboseBFDPairs';
%plist_fname = '/benchmark/PairList/Verbose3DREFPairs';

fid = fopen(plist_fname);
if fid == -1
    error('?');
end
plist = textscan(fid,'%s %s');
if size(plist,1) ~= 1 || size(plist,2) ~= 2
    error('?');
end

p1 = plist{1}; p2 = plist{2};


% compute errors i --> 1
fprintf('Computing the errors for i--> 1\n');
fprintf('0--pairwise; 1--groupwise; 2--pair icp; 3--group icp \n');

isSaveMapsMode = 0;

for i=1:Nsets
    for j=1:Nsets
        if i~=j
            name1 = dataList{j};
            name2 = dataList{i};
            
            l1 = ismember(p1, name1);
            l2 = ismember(p2, name2);
            
            
          
            if isSaveMapsMode == 1
                filename = ['benchmark/FaustMaps/', name1, '_to_',name2,'.mat'];
                
                Fmaps{i,j}{1,1}.mesh1 = loadMeshLB( name1, Ndims);
                Fmaps{i,j}{1,1}.mesh2 = loadMeshLB( name2, Ndims);
                
                map = Fmaps{i,j}{1,1}.maps{1};
                toTest = [1:Fmaps{i,j}{1,1}.mesh1.nv];
                basis1 = Fmaps{i,j}{1,1}.basis1;
                basis2 = Fmaps{i,j}{1,1}.basis2;
                dest = mapClosestDelta(map, toTest, basis1, basis2); 
                corres = [toTest; dest];
                
                save( filename, 'corres');
                
                            
            elseif sum(l1) > 0 && sum(l2) >0 && sum(l1==l2) > 0
                
                Fmaps{i,j}{1,1}.mesh1 = loadMeshLB( name1, Ndims);
                Fmaps{i,j}{1,1}.mesh2 = loadMeshLB( name2, Ndims);
                
                fprintf('Error of Fmaps %s to %s\n', name1, name2);
                % Note: Fmaps{i,j} is the mapping from j to i, so mesh1 == j,
                % mesh2 == i.
                
                nBoth = min([Fmaps{i,j}{1,1}.mesh1.nv Fmaps{i,j}{1,1}.mesh2.nv]);
                Fmaps{i,j} = mapErrorClosestDelta(Fmaps{i,j}{1,1}.mesh1, Fmaps{i,j}{1,1}.basis1,...
                    Fmaps{i,j}{1,1}.mesh2, Fmaps{i,j}{1,1}.basis2, Fmaps{i,j},1:nBoth,...
                    nn_method);
                
                if isShow == 1
                    visualizeMaps2(Fmaps{i,j}{1,1}.mesh1, Fmaps{i,j}{1,1}.basis1, Fmaps{i,j}{1,1}.mesh2,...
                        Fmaps{i,j}{1,1}.basis2, Fmaps{i,j});
                end
                
                % show figures compact
                if isShow == 2
                    
                    compactVisualizeMaps2(Fmaps{i,j}{1,1}.mesh1, Fmaps{i,j}{1,1}.basis1, Fmaps{i,j}{1,1}.mesh2,...
                        Fmaps{i,j}{1,1}.basis2, Fmaps{i,j});
                end
                
                printMapsErrors(Fmaps{i,j}{1,1}.options, Fmaps{i,j}{1,1}.mesh1, Fmaps{i,j}{1,1}.mesh2, Fmaps{i,j});
                
                % save results
                if toSave == 1
                    destinations = Fmaps{i,j}{1}.destinations;
                    fmaps = Fmaps{i,j}{1}.maps;
                    errors = Fmaps{i,j}{1}.errors;
                    %geoDistNorm = Fmaps{i,j}{1}.geoDistNormalized;
                    options = Fmaps{i,j}{1}.options;
                    numEig = Fmaps{i,j}{1}.numEig;
                    segMatches = Fmaps{i,j}{1}.segMatches;
                    segMethod = Fmaps{i,j}{1}.segMethod;
                    gwICPs = Ticp;
                    %dlmwrite(p2p_fname, (destinations-1)');
                    fmap_fname = sprintf('benchmark/BFDFmaps/eval_landmark_only/%s_to_%sN10.mat', Fmaps{i,j}{1}.mesh1.name, Fmaps{i,j}{1}.mesh2.name);
                    save(fmap_fname, 'fmaps','destinations','segMatches','options','numEig','segMethod','errors','gwICPs');
                end
                
                
            end
        end
    end
    
    out = Fmaps;
    
    
end