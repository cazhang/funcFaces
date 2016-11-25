clc;
close all;
clear all;

addpath(genpath('/local/sda7/cz679/myToolbox/manopt'));

% Purpose:
%groupwise func. map implementation proposed by W. Smith, C. Zhang

% This file is based on Manopt: www.manopt.org.
% Original author: Chao Zhang, June 4, 2015.
% Version 1.1: Chao Zhang, Aug. 3, 2015.
% Contributors: Arnaud
% Change log: evaluation on TOSCA, using geodesic error, produce CMC curve
% This is for PC demo, not for server running.
Nclass = 9;

useData = 1; % 1: cat
% set the number of shapes in each class
NclassSets = [11 6 7 9 4 8 20 12 3];
% load file name list
fileID = fopen('TOSCA_list.txt');
nameList = textscan( fileID, '%s');
fclose(fileID);

% set the basis dimension
Ndims = 30;
% Init solution
% ICP iterations
nICP = 5;

for classID = 1:Nclass
    Nsets = NclassSets(classID);
    %Nsets = 3;
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
            if i~=j && i < j
                name1 = dataList{i};
                name2 = dataList{j};
                Fmaps{i,j} = getMapofTwoShapes( name2, name1, Ndims, useData );
            end
        end
    end
    
    
    [Ndims, Nobs] = size(Fmaps{1,2}{1,1}.F);
    Npts = size(Fmaps{1,2}{1,1}.basis1, 1);
    
    % Create the problem structure
    M = stiefelfactory(Ndims, Ndims, Nsets-1); % Product of Stiefel manifold
    M_euc = euclideanfactory(Ndims, Ndims); % Euclidean manifold
    Mn = powermanifold(M_euc, Nsets-1); % product of Euc. manifold
    % P is descriptor constraints
    P   = zeros(Ndims, Nobs, Nsets);
    % Q is operator constraints
    Q   = zeros(Ndims, Ndims, Nsets);
    % T should be the mapping matrix
    T   = zeros(Ndims, Ndims, Nsets, Nsets);
    % V is the laplace beltrame basis
    V   = zeros(Ndims, Npts, Nsets); % basis
    
    for i = 1:Nsets
        if i==1
            P(:,:,i) = Fmaps{1,2}{1,1}.G;
            Q(:,:,i) = Fmaps{1,2}{1,1}.R;
            V(:,:,i) = Fmaps{1,2}{1,1}.basis2';
            
        else
            P(:,:,i) = Fmaps{1,i}{1,1}.F;
            Q(:,:,i) = Fmaps{1,i}{1,1}.S;
            V(:,:,i) = Fmaps{1,i}{1,1}.basis1';
        end
    end
    
    % set groupwise func. map problem
    gwFunmapProblem.M = M;
    % use pairwise func.map as init
    X0 = gwFunmapProblem.M.rand();
    for i=1:Nsets-1
        X0(:,:,i) = Fmaps{1,i+1}{1,1}.maps{1,1};
        %X0(:,:,i) = eye(Ndims);
    end
    % get full mapping matrix from manifold variable
    T0 = M2T( X0 );
    
    % set funmap problem
    gwFunmapProblem.cost = @(X) gwFunmapCost(X, P, Q, Nsets);
    gwFunmapProblem.grad = @(X) gwFunmapProblem.M.egrad2rgrad(X, gwFunmapGrad(X, P, Q, Nsets, Ndims));
    
    
    % Pick an algorithm to solve the problem
    options_GWFunc = [];
    options_GWFunc.maxiter = 100;
    %options.tolgradnorm = 1e-20;
    
    % groupwise func. map optimization
    % minimize ()
    
    t = tic;
    [Xsol, costopt, info] = trustregions(gwFunmapProblem, X0, options_GWFunc);
    fprintf('%g\n', toc(t));
    
    
    %% reorder basis
    %Vout = mapClosestDeltaMutu(Xsol, V, Nsets, 0);
    %[NNbasis, NNidx] = mapClosestDeltaAll(Xsol, V, Nsets);
    %figure; plot(NNidx(:,1,2));
    Vout = V;
    options_GWicp = [];
    options_GWicp.maxiter = 100;
    gwIcpProblem.M = M;
    gwIcpProblem.cost = @(X) gwIcpCost(X, Vout, Nsets, NNbasis);
    gwIcpProblem.grad = @(X) gwIcpProblem.M.egrad2rgrad(X, gwIcpGrad(X, Vout, Nsets, Ndims, NNbasis));
    
    
    %Xin = Xsol;
    Xold = Xsol;
    %Xin = initEye(Ndims, Nsets);
    %Xin = gwIcpProblem.M.rand();
    if nICP > 0
        Xicp = cell(nICP, 1);
        for i=1:nICP
            fprintf(1,'\n Group ICP Step %d \n',i);
             % nn search
            [NNbasis, NNidx] = mapClosestDeltaAll(Xold, V, Nsets);
            %figure; plot(NNidx(:,1,2));
            gwIcpProblem.cost = @(X) gwIcpCost(X, Vout, Nsets, NNbasis);
            gwIcpProblem.grad = @(X) gwIcpProblem.M.egrad2rgrad(X, gwIcpGrad(X, Vout, Nsets, Ndims, NNbasis));
            
            [Xout, costopt, info] = trustregions(gwIcpProblem, Xold, options_GWicp);
            % save
            Xicp{i} = Xout;
            Xold = Xout;
           
        end
    else
        Xicp = cell(1,1);
        Xicp{1} = Xold;
    end
    
    % display results
    isShowComparison = 2;
    toSave = 1;
    [FmapsOut] = displayResult( Fmaps, X0, Xsol, Xicp, Nsets, Ndims, nICP, isShowComparison, toSave);
    
    clear Fmaps
    clear FmapsOut
end
