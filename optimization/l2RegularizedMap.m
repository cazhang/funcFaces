function matrix = l2RegularizedMap(F, G, lambda)
% Computes a map matrix by solving:
%   min ||AF-G||_2^2 + lambda^2 * ||A||_fro^2

% Again, AF = G <--> F'A' = G'.  We augment with lambda * identity to get
% Tikhonov regularization.

lhs = [ F'; lambda * speye(size(F',2)) ];
rhs = [ G' ; sparse(size(F',2),size(G',2))];

matrix = (lhs\rhs)';