function [ssmp] = sum(x, dim)

% [S] = SUM(X, DIM) computes the sum, across all cells in x along 
% the dimension dim.
% 
% X should be an linear cell-array of matrices for which the size in at 
% least one of the dimensions should be the same for all cells 

D = checkinput(x);
if nargin==1,
  dim = find(D,1,'last');
elseif ischar(dim) && strcmp(dim, 'cell')
  % sum the cells
  if all(D)
    ssmp = zeros(size(x{1}));
    for k = 1:numel(x)
      ssmp = ssmp+x{k};
    end
  else
    error('cell-wise summing is not possible due to different dimensionality in the individual cells');
  end
  return;
else
  % check whether the requested dim can be used
  if ~D(dim), error('data can not be concatenated across dimension %d',dim); end
end
nx   = max(size(x));
ssmp = cellfun(@sum,   x, repmat({dim},1,nx), 'UniformOutput', 0);
ssmp = sum(cell2mat(ssmp), dim);  

