function [z, sd, m] = cellzscore(x, dim, flag)

% [Z, SD] = CELLZSCORE(X, DIM, FLAG) computes the zscore, across all cells in x along 
% the dimension dim, normalising by the total number of samples 
% 
% X should be an linear cell-array of matrices for which the size in at 
% least one of the dimensions should be the same for all cells. If flag==1, the mean will
% be subtracted first (default behaviour, but to save time on already demeaned data, it
% can be set to 0). SD is a vector containing the standard deviations, used for the normalisation.

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellstd');
end

if nargin<2,
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; %let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   error('no dimension to compute mean for');
  end
elseif nargin==2,
  flag = 1;
end

if flag,
  m    = cellmean(x, dim);
  x    = cellvecadd(x, -m);
end

sd   = cellstd(x, dim, 0);
z    = cellvecmult(x, 1./sd);
