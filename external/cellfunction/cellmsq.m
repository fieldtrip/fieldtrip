function [m] = cellmsq(x, dim)

% [M] = CELLMSQ(X, DIM) computes the mean of squares, across all cells in x along 
% the dimension dim.
% 
% X should be an linear cell-array of matrices for which the size in at 
% least one of the dimensions should be the same for all cells 

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellmsq');
end

if nargin==1,
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; %let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   error('no dimension to compute mean for');
  end
end

nx   = max(nx);
nsmp = cellfun('size', x, dim);
ssmp = cellfun(@sumsq,   x, repmat({dim},1,nx), 'UniformOutput', 0);
m    = sum(cell2mat(ssmp), dim)./sum(nsmp);  

function [s] = sumsq(x, dim)

s = sum(x.^2, dim);
