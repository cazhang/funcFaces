% Demo: groupwise func. map implementation proposed by C. Zhang, and W. Smith
% Reference: Zhang, Chao, et al. "Functional Faces: Groupwise Dense
% Correspondence using Functional Maps." CVPR 2016
% Original author: Chao Zhang, June 4, 2015.
% Version 1.1: C. Zhang, Aug. 3, 2015.
% Version 1.2: C. Zhang, June. 1, 2016
% Change log: this is for demo purpose

clc;close all;clear all;

path(path,'util/');
path(path,'optimization/');
path(path,'descriptors/');
path(path,'operators/');
path(path,'shape/'); 
path(path,'groupwise/');
path(path,'external/lb_code/');
path(path,'external/lb_code/MeshLP');
path(path,'external/matlab_bgl');
path(path,'external/ann_mwrapper');
path(path,'external/voronoi-laplace');
path(path,'external/voronoi-laplace/functions');
path(path,'external/voronoi-laplace/GLtree3DMex');
path(path,'external');
path(path,'cache/');
addpath(genpath('/local/sda7/cz679/myToolbox/manopt'));
%setup your path to manopt

% load file name list
fileID = fopen('data_list.txt');
nameList = textscan( fileID, '%s');
fclose(fileID);
% set nn search method
nn_method = 'ann';
% set the basis dimension
Ndims = 30;
% set iters for ICP
nICP = 10;
classID = 1;
% number of data
Nsets = 3;
Fmaps = cell(Nsets, Nsets);
dataList = cell(Nsets,1);

% load shape list
for objID = 1:Nsets
    if classID == 1
        dataList{objID} = nameList{1}{objID};
    else
        startID = sum(NclassSets(1:classID-1));
        dataList{objID} = nameList{1}{startID + objID};
    end
end
% init with pariwise func. map
for i=1:Nsets
    for j=1:Nsets
        if i==j
            Fmaps{i,j}{1,1}.maps = eye(Ndims);
        elseif i~=j
            name1 = dataList{i};
            name2 = dataList{j};
            fmap_name = [name2,'_to_',name1,'_d',num2str(Ndims),'_10pt.mat'];
            if ~exist(fmap_name, 'file')
                % compute pairwise map
                Fmaps{i,j} = getMapofTwoFaceShapes( name2, name1, Ndims );
                maps = Fmaps{i,j}{1,1}.maps;
                ppMaps = Fmaps{i,j}{1,1}.ppMaps;
                save( fmap_name, 'maps', 'ppMaps');
            else
                fprintf('loading %s->%s\n', name2, name1 );
                % load pairwise map
                Fmaps{i,j} = fillMapofTwoFaceShapes( name2, name1, Ndims );
                load(fmap_name);
                Fmaps{i,j}{1,1}.maps = maps;
                Fmaps{i,j}{1,1}.ppMaps = ppMaps;
                
            end
        end
    end
end

[Ndims, Nobs] = size(Fmaps{1,2}{1,1}.F);
Npts = size(Fmaps{1,2}{1,1}.basis1, 1);

% Create the problem structure: Product of Stiefel manifold
M = stiefelfactory(Ndims, Ndims, Nsets-1);
% P is descriptor constraints
P   = zeros(Ndims, Nobs, Nsets);
% Q is operator constraints
Q   = zeros(Ndims, Ndims, Nsets);
% T should be the mapping matrix
T   = zeros(Ndims, Ndims, Nsets, Nsets);
% V is the LB basis
V   = zeros(Ndims, Npts, Nsets); % basis

V = cell(1,Nsets);
for i = 1:Nsets
    if i==1
        P(:,:,i) = Fmaps{1,2}{1,1}.G;
        Q(:,:,i) = Fmaps{1,2}{1,1}.R;
        V{1,i} = Fmaps{1,2}{1,1}.basis2';
        
    else
        P(:,:,i) = Fmaps{1,i}{1,1}.F;
        Q(:,:,i) = Fmaps{1,i}{1,1}.S;
        V{1,i} = Fmaps{1,i}{1,1}.basis1';
    end
end

% set groupwise func. map problem
gwFunmapProblem.M = M;
% use pairwise func.map as init
X0 = gwFunmapProblem.M.rand();
X1 = gwFunmapProblem.M.rand();

for i=1:Nsets-1
    % X0: pairwise func map
    X0(:,:,i) = Fmaps{1,i+1}{1,1}.maps{1,1};
    % X1: pairwise+icp func map
    X1(:,:,i) = Fmaps{1,i+1}{1,1}.ppMaps{1,1};
end

% set funmap problem
gwFunmapProblem.cost = @(X) gwFunmapCost(X, P, Q, Nsets);
gwFunmapProblem.grad = @(X) gwFunmapProblem.M.egrad2rgrad(X, gwFunmapGrad(X, P, Q, Nsets, Ndims));

% Pick an algorithm to solve the problem
options_GWFunc = [];
options_GWFunc.maxiter = 60;
%options.tolgradnorm = 1e-20;
Xsol = cell(3,1);
% pairwise func.map init
fprintf('Starting Groupwise Functional Map Optimisation...')
t = tic;
[Xsol{1}, costopt, info] = trustregions(gwFunmapProblem, X0, options_GWFunc);
fprintf('Group Running time is: %g\n', toc(t));

%% Group ICP process
options_GWicp = [];
options_GWicp.maxiter = 30;
options_GWicp.step = 10;
gwIcpProblem.M = M;

Xold1 = Xsol{1};

fprintf('Starting Groupwise-ICP Functional Map Optimisation...')
if nICP > 0
    Xicp = cell(nICP, 1);
    for i=1:nICP
        fprintf(1,'\n Group ICP Step %d \n',i);
        
        %% single direction NN
        [NNbasis1] = mapClosestDeltaAllHeter(Xold1, V, Nsets, options_GWicp.step);
        
        gwIcpProblem.cost = @(X) gwIcpCostHeter(X, Nsets, NNbasis1);
        gwIcpProblem.grad = @(X) gwIcpProblem.M.egrad2rgrad(X, gwIcpGradHeter(X, Nsets, Ndims, NNbasis1));

        [Xold1, costopt, info] = trustregions(gwIcpProblem, Xold1, options_GWicp);
        Xicp{i} = Xold1;
    end
end


%% compute NN indices
nnidx_map = cell(Nsets, Nsets);
nnidx_icp = cell(Nsets, Nsets);
nnidx_gru = cell(Nsets, Nsets);
nnidx_gru2 = cell(Nsets, Nsets);

for j=1:Nsets
    for i=1:Nsets
        if i~=j
            % j -> i
            fprintf('i = %d, j = %d\n', i, j);

            Ticp = M2T(X1); map = Ticp(:,:,i,j);
            pts1 = map * Fmaps{i,j}{1,1}.basis1';
            pts2 = Fmaps{i,j}{1,1}.basis2';
            t = tic;
            nnidx_icp{i,j} = annquery(pts2, pts1, 1);
            fprintf(' ICP %g\n', toc(t));
            
            Tgru2 = M2T(Xicp{nICP}); map = Tgru2(:,:,i,j);
            pts1 = map * Fmaps{i,j}{1,1}.basis1';
            pts2 = Fmaps{i,j}{1,1}.basis2';
            t = tic;
            nnidx_gru2{i,j} = annquery(pts2, pts1, 1);
            fprintf('Group ICP %g\n', toc(t));
           
        end
    end
end


%% show ICP results
j = 1;
i = 2;
k = 3;
h = figure;
mesh1 = loadMeshLB(dataList{j}, Ndims);
mesh2 = loadMeshLB(dataList{i}, Ndims);
mesh3 = loadMeshLB(dataList{k}, Ndims);

X1 = mesh1.vertices; T1 = mesh1.triangles;
X2 = mesh2.vertices; T2 = mesh2.triangles;
X3 = mesh3.vertices; T3 = mesh3.triangles;

%%
col1 = mesh1.texture;

h1 = subplot(331);
patch('Faces',T1,'Vertices',X1,'FaceColor','interp', ...
    'FaceVertexCData', col1, 'EdgeColor', 'none');
axis equal; axis tight; axis off; cameratoolbar;

h2 = subplot(332);
col2 = mesh2.texture;
patch('Faces',T2,'Vertices',X2,'FaceColor','interp', ...
    'FaceVertexCData', col2, 'EdgeColor', 'none');
axis equal; axis tight; axis off; cameratoolbar;

h3 = subplot(333);
col3 = mesh3.texture;
patch('Faces',T3,'Vertices',X3,'FaceColor','interp', ...
    'FaceVertexCData', col3, 'EdgeColor', 'none');
axis equal; axis tight; axis off; cameratoolbar;

h4 = subplot(334);
col4 = mesh2.texture(nnidx_icp{i,j},:);
patch('Faces',T1,'Vertices',X1,'FaceColor','interp', ...
    'FaceVertexCData', col4, 'EdgeColor', 'none');
axis equal; axis tight; axis off; cameratoolbar;

h5 = subplot(335);
col5 = mesh3.texture(nnidx_icp{k,j},:);
patch('Faces',T1,'Vertices',X1,'FaceColor','interp', ...
    'FaceVertexCData', col5, 'EdgeColor', 'none');
axis equal; axis tight; axis off; cameratoolbar;

col6 = var([col1(:,1) col4(:,1) col5(:,1)]');
col6 = col6'./max(col6);

col7 = var([col1(:,2) col4(:,2) col5(:,2)]');
col7 = col7'./max(col7);

col8 = var([col1(:,3) col4(:,3) col5(:,3)]');
col8 = col8'./max(col8);

col9 = [col6 col7 col8];
col10 = mean(col9, 2);
fprintf('icp variance = %f\n', mean(col10));
h6 = subplot(336);

patch('Faces',T1,'Vertices',X1,'FaceColor','interp', ...
    'FaceVertexCData', col10, 'EdgeColor', 'none');
axis equal; axis tight; axis off; cameratoolbar;
colormap(jet(256));

h7= subplot(337);

col11 = mesh2.texture(nnidx_gru2{i,j},:);
patch('Faces',T1,'Vertices',X1,'FaceColor','interp', ...
    'FaceVertexCData', col11, 'EdgeColor', 'none');
axis equal; axis tight; axis off; cameratoolbar;

h8 = subplot(338);

col12 = mesh3.texture(nnidx_gru2{k,j},:);
patch('Faces',T1,'Vertices',X1,'FaceColor','interp', ...
    'FaceVertexCData', col12, 'EdgeColor', 'none');
axis equal; axis tight; axis off; cameratoolbar;

col13 = var([col1(:,1) col11(:,1) col12(:,1)]');
col13 = col13'./max(col13);

col14 = var([col1(:,2) col11(:,2) col12(:,2)]');
col14 = col14'./max(col14);

col15 = var([col1(:,3) col11(:,3) col12(:,3)]');
col15 = col15'./max(col15);

col16 = [col13 col14 col15];
col17 = mean(col16, 2);
fprintf('gicp variance = %f\n', mean(col17));

h9 = subplot(339);
patch('Faces',T1,'Vertices',X1,'FaceColor','interp', ...
    'FaceVertexCData', col17, 'EdgeColor', 'none');
axis equal; axis tight; axis off; cameratoolbar;
colormap(jet(256));
set(h,'WindowStyle','normal');
