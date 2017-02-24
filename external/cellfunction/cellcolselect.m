function [y] = cellcolselect(x, cols)

% [Y] = CELLCOLSELECT(X, COLS) outputs cell-array Y with the same dimensionality as X
% Each cell in Y only contains the columns COLS from the original corresponding cell in X
% 
% X (and Y) should be linear cell-array(s) of matrices for which the number of rows 
% should be the same for all cells

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellrows');
end

y = cellfun(@colc, x, repmat(mat2cell(cols(:),length(cols),1),nx), 'UniformOutput', 0);

function [y] = colc(x, cols)

y = x(:,cols);
