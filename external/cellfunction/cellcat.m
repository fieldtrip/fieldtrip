function [z] = cellcat(dim, x, y)

% [Z] = CELLCAT(DIM, X, Y) outputs cell-array Z with the same dimensionality as X and Y
% Each cell in Z contains a concatenation of the corresponding arrays in X and Y,
% along dimension DIM
% 
% X (and Y) should be linear cell-array(s) of matrices for which the number of rows (if dim=1)
% or columns (if dim=2) should be the same for all cells

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellcat');
end

z = cellfun(@catc, x, y, repmat({dim},nx), 'UniformOutput', 0);

function [z] = catc(x, y, dim)

z = cat(dim, x, y);