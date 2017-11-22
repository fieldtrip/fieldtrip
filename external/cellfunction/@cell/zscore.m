function [z, mu, s] = zscore(x, varargin)

% [Z, MU, S] = ZSCORE(X, NORMALIZEFLAG, DIM, FLAG) computes the zscore, across all cells in x along 
% the dimension dim, normalising by the total number of samples 
% 
% X should be an linear cell-array of matrices for which the size in at 
% least one of the dimensions should be the same for all cells. If flag==1, the mean will
% be subtracted first (default behaviour, but to save time on already demeaned data, it
% can be set to 0). MU and S are vectors containing the mean and standard deviations used
% for zscoring

if numel(varargin)==0, varargin{1} = []; end
if numel(varargin)==1, varargin{2} = []; end
if numel(varargin)==2, varargin{3} = 1;  end

if varargin{3},
  mu = nanmean(x, varargin{2});
  x  = cellvecadd(x, -mu);
end
s = nanstd(x, varargin{1:end});
z = cellvecmult(x, 1./s);

if ~varargin{3},
  mu = zeros(size(s))+nan;
end