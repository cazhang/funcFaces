function descriptors = coeffsToDescriptors(coeffs, basis)
% This method is the reverse of descriptorsToCoeffs.

descriptors = basis * coeffs;