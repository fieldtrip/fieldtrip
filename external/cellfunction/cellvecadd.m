function [y] = cellvecadd(x, v)

% [Y]= CELLVECADD(X, V) - add vector to all rows or columns of each matrix 
% in cell-array X

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellvecadd');
end

if ~iscell(v),
  v = repmat({v}, nx);
end

y  = cellfun(@bsxfun, repmat({@plus}, nx), x, v, 'UniformOutput', 0);

