function [y] = cellrowselect(x, rows)

% [Y] = CELLROWSELECT(X, ROWS) outputs cell-array Y with the same dimensionality as X
% Each cell in Y only contains the rows ROWS from the original corresponding cell in X
% 
% X (and Y) should be linear cell-array(s) of matrices for which the number of rows 
% should be the same for all cells

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellrows');
end

y = cellfun(@rowc, x, repmat(mat2cell(rows(:),length(rows),1),nx), 'UniformOutput', 0);

function [y] = rowc(x, rows)

y = x(rows,:);
