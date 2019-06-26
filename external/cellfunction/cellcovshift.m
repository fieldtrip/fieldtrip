function [c] = cellcovshift(x, shift, dim, flag)

% [C] = CELLCOVSHIFT(X, SHIFT, DIM) computes the covariance, across all cells
% in x along the dimension dim. 
% 
% X should be linear cell-array(s) of matrices for which the size in at 
% least one of the dimensions is be the same for all cells 

if nargin<4,
  flag = 1;
end

if nargin<3,
  dim = find(size(x{1})>1, 1, 'first');
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1) || ndims(x{1})>2,
  error('incorrect input for cellcovshift');
end

% shift the time axis
y = cellshift(x, shift, dim);

% compute covariance
c = cov(y, 1, dim, flag);
%c = cellnancov(y, 1, dim, flag);