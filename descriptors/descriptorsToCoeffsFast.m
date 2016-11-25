function coeffs = descriptorsToCoeffsFast(descriptors, basisi)
% This method takes a matrix where each column is a different descriptor.
% It returns another matrix whose columns are the coefficients of the
% descriptors projected onto the provided basis.

coeffs = basisi * descriptors;
