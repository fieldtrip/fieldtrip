function [y] = cellvecmult(x, v)

% [Y]= CELLVECMULT(X, V) - multiply vectors in cell-array V
% to all rows or columns of each matrix in cell-array X
% V can be a vector or a cell-array of vectors

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellvecmult');
end

if ~iscell(v),
  v = repmat({v}, nx);
end

sx1 = cellfun('size', x, 1);
sx2 = cellfun('size', x, 2);
sv1 = cellfun('size', v, 1);
sv2 = cellfun('size', v, 2);
if all(sx1==sv1) && all(sv2==1),    
elseif all(sx2==sv2) && all(sv1==1),
elseif all(sv1==1) && all(sv2==1),
else   error('inconsistent input');
end  

y  = cellfun(@bsxfun, repmat({@times}, nx), x, v, 'UniformOutput', 0);
