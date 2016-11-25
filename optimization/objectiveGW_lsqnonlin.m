function [F, G] = objectiveGW_lsqnonlin(X, P, Q)
% P_i is the descriptor constraint, Ndims x Nobs
% Q_i is the operator constraint, Ndims x Ndims
% T_i is the mapping from i-th shape to 1st shape, Ndims x Ndims
% U_i is the mapping from 1st shape to i-th shape, Ndism x Ndims
% the objective is: sum_{ij} 


Ndims   = size(P, 1);
Nobs    = size(P, 2);
Nsets   = size(P, 3);

alpha = 1.0;

T       = cell(Nsets, 1);
U       = cell(Nsets, 1);
T{1}    = eye(Ndims);
U{1}    = eye(Ndims);
I = eye(Ndims);
for i = 2: Nsets
    % back
    T{i} = reshape(X(1 + (i - 2) * Ndims^2: (i - 1) * Ndims^2), Ndims, Ndims);
    % forward
    U{i} = inv(T{i});
    
end


F1 = zeros(Ndims, Nobs, Nsets, Nsets);

F2 = zeros(Ndims, Ndims, Nsets, Nsets);

F3 = zeros(Ndims, Ndims, Nsets, Nsets);



if nargout > 1
    G = zeros(Ndims, Ndims, Nsets - 1);
end

for i = 1: Nsets
    for j = 1: Nsets
        if i ~= j
            map = U{i} * T{j};
            % descriptor constraints
            F1(:, :, i, j) = P(:, :, i) - map * P(:, :, j);
            % operator constraints
            temp2 = Q(:, :, j) * map - map * Q(:, :, i);
            F2(:, :, i, j) = temp2 / norm(temp2, inf);
            % orthogonality
            temp3 = map' * map - I;
            F3(:, :, i, j) = alpha * temp3 / norm(temp3, inf);
            
            F = [F1(:); F2(:); F3(:)];
            
            %temp = [temp1(:); temp2(:); temp3(:)];
            %O = O + norm(temp, 2)^2;
            %O = O + norm(vec(P(:, :, i) - (U{i}(1: Ndims, :) * T{j}) * [P(:, :, j); e]), 2)^2;
         
            if nargout > 1 && i ~= 1
                tmp             = (U{i} * T{j}) * P(:, :, j);
                G(:, :, i - 1)  = G(:, :, i - 1) + U{i}.' * (P(:, :, i) - tmp) * tmp.' - U{j}.' * (P(:, :, j) - (U{j} * T{i}) * P(:, :, i)) * P(:, :, i).';
            end
        end
    end
end

clear temp temp1 temp2 temp3 map

if nargout > 1
    temp = G(:, :, :);
    temp2 = temp(:);
    G = 2 * temp2;
    %G = 2 * vec(G(1: Ndims, :, :));
    
end