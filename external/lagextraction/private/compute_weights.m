function weights = compute_weights(distances,sigma)
% compute the weights matrix from the distance matrix
% with a gaussian kernel

% $Id: compute_weights.m 2 2009-06-16 19:24:10Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-06-16 15:24:10 -0400 (Mar, 16 jui 2009) $
% $Revision: 2 $

weights = exp(-distances.^2/sigma^2);
