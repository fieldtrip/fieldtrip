function [m] = nanmean(x, dim)

% [M] = NANMEAN(X, DIM) computes the mean, across all cells in x along 
% the dimension dim, accounting for NaNs.
% 
% X should be an linear cell-array of matrices for which the size in at 
% least one of the dimensions should be the same for all cells 

% D = checkinput(x);
% if nargin==1,
%   dim = find(D,1,'last');
% else
%   % check whether the requested dim can be used
%   if ~D(dim), error('data can not be concatenated across dimension %d',dim); end
% end
if nargin < 2
  dim = [];
end
nx   = max(size(x));
nsmp = cellfun(@sum,    isfinite(x), repmat({dim},1,nx), 'UniformOutput', 0);
ssmp = cellfun(@nansum,           x, repmat({dim},1,nx), 'UniformOutput', 0);
m    = nansum(cell2mat(ssmp), dim)./sum(nsmp);

