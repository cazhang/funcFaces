function coeffs = descriptorsToCoeffs(descriptors, basis)
% This method takes a matrix where each column is a different descriptor.
% It returns another matrix whose columns are the coefficients of the
% descriptors projected onto the provided basis.

% TODO - if this is too slow, prefactor basis, or precompute pinv(basis)
coeffs = basis \ descriptors;
