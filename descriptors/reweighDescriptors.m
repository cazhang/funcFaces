function [FRef, GRef] = reweighDescriptors(F, G, nbEigen, options)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
normalization = sqrt(sum(F.^2));
FRef = F./repmat(normalization, [nbEigen, 1]);

normalization = sqrt(sum(G.^2));
GRef = G./repmat(normalization, [nbEigen, 1]);


% FRef(:, end-2:end) = FRef(:, end-2:end)*100;
% GRef(:, end-2:end) = GRef(:, end-2:end)*100;
% list = [21 22 57 20 42 72 56 19];
% list = list - 3;
% 
% for i=1:length(list)
%     start = (list(i)-1)*10+1;
%     ending = list(i)*10;
%     
%     FRef(:, start:ending) = FRef(:, start:ending)*100;
%     GRef(:, start:ending) = GRef(:, start:ending)*100;
%     
% end


