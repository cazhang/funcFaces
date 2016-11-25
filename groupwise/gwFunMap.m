%% Reset

clear all
close all
clc

%addpath('./src/');
%addpath(genpath('/local/sda7/cz679/myToolbox/yalmip'));
addpath(genpath('/local/sda7/cz679/funMap/fmaps_siggraph2012_code_updated/fmaps_siggraph2012_code'));

%% Prepare data
% descriptors for P(1)
name1 = 'cat0';
name2 = 'cat1';
name3 = 'cat5';
name4 = 'cat6';
name5 = 'cat10';

numEigen = 20;
  


%Fmaps12 = getMapofTwoShapes( name1, name2, numEigen );
%Fmaps13 = getMapofTwoShapes( name1, name3, numEigen );
%Fmaps14 = getMapofTwoShapes( name1, name4, numEigen );
Fmaps21 = getMapofTwoShapes( name2, name1, numEigen );
Fmaps31 = getMapofTwoShapes( name3, name1, numEigen );
Fmaps41 = getMapofTwoShapes( name4, name1, numEigen );
Fmaps51 = getMapofTwoShapes( name5, name1, numEigen );



%% groupwise opti
[Ndims, Nobs] = size(Fmaps21{1,1}.F);
Nsets = 3;
P   = zeros(Ndims, Nobs, Nsets);
Q   = zeros(Ndims, Ndims, Nsets);
T   = zeros(Ndims, Ndims, Nsets, Nsets);

if Nsets >=1
P(:,:,1) = Fmaps21{1,1}.G;
Q(:,:,1) = Fmaps21{1,1}.R;
T(:, :, 1, 1) = eye(Ndims);
end
if Nsets >=2
P(:,:,2) = Fmaps21{1,1}.F;
Q(:,:,2) = Fmaps21{1,1}.S;
%T(:, :, 1, 2) = Fmaps21{1,1}.maps{1,1};
%T(:, :, 1, 2) = Fmaps21{1,1}.maps{1,1};
T(:, :, 1, 2) = eye(Ndims);
end
if Nsets >=3
P(:,:,3) = Fmaps31{1,1}.F;
Q(:,:,3) = Fmaps31{1,1}.S;
%T(:, :, 1, 3) = Fmaps31{1,1}.maps{1,1};
%T(:, :, 1, 3) = Fmaps31{1,1}.maps{1,1};
T(:, :, 1, 3) = eye(Ndims);
end
if Nsets >=4
P(:,:,4) = Fmaps41{1,1}.F;
Q(:,:,4) = Fmaps41{1,1}.S;
T(:, :, 1, 4) = Fmaps41{1,1}.maps{1,1};
%T(:, :, 1, 4) = eye(Ndims);
end
if Nsets >=5
P(:,:,5) = Fmaps51{1,1}.G;
Q(:,:,5) = Fmaps51{1,1}.R;
T(:, :, 1, 5) = Fmaps51{1,1}.maps{1,1};
%T(:, :, 1, 5) = eye(Ndims);
end


%P = P + 0.01 * randn(Ndims, Nobs, Nsets);


for i = 2: Nsets
    T(:, :, i, i)   = eye(Ndims);
    %T(:, :, 1, i)   = T(:, :, i, 1)^-1;
    T(:, :, i, 1) = T(:, :, 1, i)^-1;
    
    for j = 2: i - 1
        T(:, :, i, j) = T(:, :, i, 1) * T(:, :, 1, j);
        T(:, :, j, i) = T(:, :, j, 1) * T(:, :, 1, i);
    end
end


%% Constrained formulation

options             = optimset;
options.Jacobian    = 'off';
options.Display     = 'iter';
options.MaxIter     = 50;
options.GradObj     = 'on';
% options.Jacobian    = 'off';
Omin                = inf;

%options = optimoptions('PlotFcns',@optimplotfval);
%matlabpool('open',4);
for i = 1: 1
    fprintf('Constrained optimization : Run %d of total %d\n', i, Nsets);
    %temp = That(:, :, 1, 2: Nsets);
    temp = T(:, :, 1, 2: Nsets);
    X0   = temp(:);
    
    tic
    %[X, O]  = lsqnonlin(@(X)objectiveGW_lsqnonlin(X, P, Q), X0, [], [], options);
    [X, O]  = fminunc(@(X)objectiveGW_soft(X, P, Q), X0, options);
    toc  
   
end

%matlabpool close;

That                = zeros(Ndims, Ndims, Nsets, Nsets);
That(:, :, 1, 1)    = eye(Ndims);
for i = 2: Nsets
    That(:, :, i, i) = eye(Ndims);
    That(:, :, 1, i) = reshape(X(1 + (i - 2) * Ndims^2: (i - 1) * Ndims^2), Ndims, Ndims);
    
    That(:, :, i, 1) = That(:, :, 1, i)^-1;
    %That(:, :, i, 1) = That(:, :, 1, i)';
    for j = 2: i - 1
        That(:, :, i, j) = That(:, :, i, 1) * That(:, :, 1, j);
        That(:, :, j, i) = That(:, :, j, 1) * That(:, :, 1, i);
    end
end

%That = circshift(That, [0, 0, 1 - imin, 1 - imin]);

clear i j k imin Omin Xmin X X0 options

%% display before
fprintf('Objective value: %.6e\n', O);

figure;
for i = 1: Nsets
    for j = 1: Nsets
        subplot(Nsets, Nsets, i + (j - 1) * Nsets)
        %tmp = (T(:, :, i, j) - That(:, :, i, j)).^2;
        tmp = T(:, :, i, j);
        imagesc(tmp), colorbar
        tmp = sum(tmp(:)) / (Ndims * Ndims);
        if tmp == 0
            caxis([0 1]);
        end
        title(sprintf('That_{%d%d} vs T_{%d%d} (%.2e)', i, j, i, j, tmp));
        fprintf('Estimation error %d-%d: %.2e\n', i, j, tmp);
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
    end
end

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

%clear i j k tmp O

%% display after
fprintf('Objective value: %.6e\n', O);

figure;
for i = 1: Nsets
    for j = 1: Nsets
        subplot(Nsets, Nsets, i + (j - 1) * Nsets)
        %tmp = (T(:, :, i, j) - That(:, :, i, j)).^2;
        tmp = That(:, :, i, j);
        imagesc(tmp), colorbar
        tmp = sum(tmp(:)) / (Ndims * Ndims);
        if tmp == 0
            caxis([0 1]);
        end
        title(sprintf('That_{%d%d} vs T_{%d%d} (%.2e)', i, j, i, j, tmp));
        fprintf('EstimationO error %d-%d: %.2e\n', i, j, tmp);
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
    end
end

for i = 1: Nsets
    for j = 1: Nsets
        if i ~= j
            tmp = (That(:, :, i, j) * That(:, :, j, i) - eye(Ndims)).^2;
            fprintf('Invertibility constraint %d-%d: %.2e\n', i, j, sum(tmp(:)) / (Ndims * Ndims))
        end
    end
end

for i = 1: Nsets
    for j = 1: Nsets
        for k = 1: Nsets
            if i ~= j && j ~= k && i ~= k
                tmp = (That(:, :, i, j) * That(:, :, j, k) - That(:, :, i, k)).^2;
                fprintf('Transitivity constraint %d-%d-%d: %.2e\n', i, j, k, sum(tmp(:)) / (Ndims * Ndims))
            end
        end
    end
end

clear i j k tmp O


%% compare pairwise maps with groupwise maps
% only compare i->1 mapping

if Nsets >=2
    Fmaps21{1,1}.maps{2} = That(:, :, 1, 2);
    Fmaps21{1,1}.params{2} = 0;
end
if Nsets >=3
    Fmaps31{1,1}.maps{2} = That(:, :, 1, 3);
    Fmaps31{1,1}.params{2} = 0;
end
if Nsets >=4
    Fmaps41{1,1}.maps{2} = That(:, :, 1, 4);
    Fmaps41{1,1}.params{2} = 0;
end
if Nsets >=5
    Fmaps51{1,1}.maps{2} = That(:, :, 1, 5);
    Fmaps51{1,1}.params{2} = 0;
end

% compute i -> 1
if Nsets >=2
    Fmaps12{1,1} = Fmaps21{1,1};
    Fmaps12{1,1}.mesh1 = Fmaps21{1,1}.mesh2;
    Fmaps12{1,1}.mesh2 = Fmaps21{1,1}.mesh1;
    Fmaps12{1,1}.basis1 = Fmaps21{1,1}.basis2;
    Fmaps12{1,1}.basis2 = Fmaps21{1,1}.basis1;
    
    Fmaps12{1,1}.maps{1} = T(:, :, 2, 1);
    Fmaps12{1,1}.params{1} = 0;
    Fmaps12{1,1}.maps{2} = That(:, :, 2, 1);
    Fmaps12{1,1}.params{2} = 0;
end
if Nsets >=3
    Fmaps13{1,1} = Fmaps31{1,1};
    Fmaps13{1,1}.mesh1 = Fmaps31{1,1}.mesh2;
    Fmaps13{1,1}.mesh2 = Fmaps31{1,1}.mesh1;
    Fmaps13{1,1}.basis1 = Fmaps31{1,1}.basis2;
    Fmaps13{1,1}.basis2 = Fmaps31{1,1}.basis1;
    
    Fmaps13{1,1}.maps{1} = T(:, :, 3, 1);
    Fmaps13{1,1}.params{1} = 0;
    Fmaps13{1,1}.maps{2} = That(:, :, 3, 1);
    Fmaps13{1,1}.params{2} = 0;
end
if Nsets >=4
    Fmaps14{1,1} = Fmaps41{1,1};
    Fmaps14{1,1}.mesh1 = Fmaps41{1,1}.mesh2;
    Fmaps14{1,1}.mesh2 = Fmaps41{1,1}.mesh1;
    Fmaps14{1,1}.basis1 = Fmaps41{1,1}.basis2;
    Fmaps14{1,1}.basis2 = Fmaps41{1,1}.basis1;
    
    Fmaps14{1,1}.maps{1} = T(:, :, 4, 1);
    Fmaps14{1,1}.params{1} = 0;
    Fmaps14{1,1}.maps{2} = That(:, :, 4, 1);
    Fmaps14{1,1}.params{2} = 0;
end
if Nsets >=5
    Fmaps15{1,1} = Fmaps51{1,1};
    Fmaps15{1,1}.mesh1 = Fmaps51{1,1}.mesh2;
    Fmaps15{1,1}.mesh2 = Fmaps51{1,1}.mesh1;
    Fmaps15{1,1}.basis1 = Fmaps51{1,1}.basis2;
    Fmaps15{1,1}.basis2 = Fmaps51{1,1}.basis1;
    
    Fmaps15{1,1}.maps{1} = T(:, :, 5, 1);
    Fmaps15{1,1}.params{1} = 0;
    Fmaps15{1,1}.maps{2} = That(:, :, 5, 1);
    Fmaps15{1,1}.params{2} = 0;
end

% compute errors i --> 1
fprintf('Computing the errors\n');
if Nsets >= 2
    fprintf('Error of Fmaps 2->1...\n');
    Fmaps21 = mapErrorClosestDelta(Fmaps21{1,1}.mesh1, Fmaps21{1,1}.basis1,...
    Fmaps21{1,1}.mesh2, Fmaps21{1,1}.basis2, Fmaps21, 1:(Fmaps21{1,1}.mesh1.nv));
end
if Nsets >= 3
    fprintf('Error of Fmaps 3->1...\n');
    Fmaps31 = mapErrorClosestDelta(Fmaps31{1,1}.mesh1, Fmaps31{1,1}.basis1,...
    Fmaps31{1,1}.mesh2, Fmaps31{1,1}.basis2, Fmaps31, 1:(Fmaps31{1,1}.mesh1.nv));
end
if Nsets >= 4
    fprintf('Error of Fmaps 1->4...\n');
    Fmaps41 = mapErrorClosestDelta(Fmaps41{1,1}.mesh1, Fmaps41{1,1}.basis1,...
    Fmaps41{1,1}.mesh2, Fmaps41{1,1}.basis2, Fmaps41, 1:(Fmaps41{1,1}.mesh1.nv));
end
if Nsets >= 5
    fprintf('Error of Fmaps 1->5...\n');
    Fmaps51 = mapErrorClosestDelta(Fmaps51{1,1}.mesh1, Fmaps51{1,1}.basis1,...
    Fmaps51{1,1}.mesh2, Fmaps51{1,1}.basis2, Fmaps51, 1:(Fmaps51{1,1}.mesh1.nv));
end

% % Visualize
if Nsets >= 2
    visualizeMaps2(Fmaps21{1,1}.mesh1, Fmaps21{1,1}.basis1, Fmaps21{1,1}.mesh2,...
    Fmaps21{1,1}.basis2, Fmaps21);
    printMapsErrors(Fmaps21{1,1}.options, Fmaps21{1,1}.mesh1, Fmaps21{1,1}.mesh2, Fmaps21);
end

if Nsets >= 3
    visualizeMaps2(Fmaps31{1,1}.mesh1, Fmaps31{1,1}.basis1, Fmaps31{1,1}.mesh2,...
    Fmaps31{1,1}.basis2, Fmaps31);
    printMapsErrors(Fmaps31{1,1}.options, Fmaps31{1,1}.mesh1, Fmaps31{1,1}.mesh2, Fmaps31);
end

if Nsets >= 4
    visualizeMaps2(Fmaps41{1,1}.mesh1, Fmaps41{1,1}.basis1, Fmaps41{1,1}.mesh2,...
    Fmaps41{1,1}.basis2, Fmaps41);
    printMapsErrors(Fmaps41{1,1}.options, Fmaps41{1,1}.mesh1, Fmaps41{1,1}.mesh2, Fmaps41);
end

if Nsets >= 5
    visualizeMaps2(Fmaps51{1,1}.mesh1, Fmaps51{1,1}.basis1, Fmaps51{1,1}.mesh2,...
    Fmaps51{1,1}.basis2, Fmaps51);
    printMapsErrors(Fmaps51{1,1}.options, Fmaps51{1,1}.mesh1, Fmaps51{1,1}.mesh2, Fmaps51);
end


% compute errors 1 --> i
fprintf('Computing the errors\n');
if Nsets >= 2
    fprintf('Error of Fmaps 2->1...\n');
    Fmaps12 = mapErrorClosestDelta(Fmaps12{1,1}.mesh1, Fmaps12{1,1}.basis1,...
    Fmaps12{1,1}.mesh2, Fmaps12{1,1}.basis2, Fmaps12, 1:(Fmaps12{1,1}.mesh1.nv));
end
if Nsets >= 3
    fprintf('Error of Fmaps 3->1...\n');
    Fmaps13 = mapErrorClosestDelta(Fmaps13{1,1}.mesh1, Fmaps13{1,1}.basis1,...
    Fmaps13{1,1}.mesh2, Fmaps13{1,1}.basis2, Fmaps13, 1:(Fmaps13{1,1}.mesh1.nv));
end
if Nsets >= 4
    fprintf('Error of Fmaps 1->4...\n');
    Fmaps14 = mapErrorClosestDelta(Fmaps14{1,1}.mesh1, Fmaps14{1,1}.basis1,...
    Fmaps14{1,1}.mesh2, Fmaps14{1,1}.basis2, Fmaps14, 1:(Fmaps14{1,1}.mesh1.nv));
end
if Nsets >= 5
    fprintf('Error of Fmaps 1->5...\n');
    Fmaps15 = mapErrorClosestDelta(Fmaps15{1,1}.mesh1, Fmaps15{1,1}.basis1,...
    Fmaps15{1,1}.mesh2, Fmaps15{1,1}.basis2, Fmaps15, 1:(Fmaps15{1,1}.mesh1.nv));
end

% % Visualize
if Nsets >= 2
    visualizeMaps2(Fmaps12{1,1}.mesh1, Fmaps12{1,1}.basis1, Fmaps12{1,1}.mesh2,...
    Fmaps12{1,1}.basis2, Fmaps12);
    printMapsErrors(Fmaps12{1,1}.options, Fmaps12{1,1}.mesh1, Fmaps12{1,1}.mesh2, Fmaps12);
end

if Nsets >= 3
    visualizeMaps2(Fmaps13{1,1}.mesh1, Fmaps13{1,1}.basis1, Fmaps13{1,1}.mesh2,...
    Fmaps13{1,1}.basis2, Fmaps13);
    printMapsErrors(Fmaps13{1,1}.options, Fmaps13{1,1}.mesh1, Fmaps13{1,1}.mesh2, Fmaps13);
end

if Nsets >= 4
    visualizeMaps2(Fmaps14{1,1}.mesh1, Fmaps14{1,1}.basis1, Fmaps14{1,1}.mesh2,...
    Fmaps14{1,1}.basis2, Fmaps14);
    printMapsErrors(Fmaps14{1,1}.options, Fmaps14{1,1}.mesh1, Fmaps14{1,1}.mesh2, Fmaps14);
end

if Nsets >= 5
    visualizeMaps2(Fmaps15{1,1}.mesh1, Fmaps15{1,1}.basis1, Fmaps15{1,1}.mesh2,...
    Fmaps15{1,1}.basis2, Fmaps15);
    printMapsErrors(Fmaps15{1,1}.options, Fmaps15{1,1}.mesh1, Fmaps15{1,1}.mesh2, Fmaps15);
end