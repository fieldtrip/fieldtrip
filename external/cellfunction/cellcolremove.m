function [y] = cellcolremove(x, cols)

% [Y] = CELLCOLREMOVE(X, COLS) outputs cell-array Y with the same dimensionality as X
% Each cell in Y has the indexed columns COLS removed from the original corresponding cell in X
% 
% X (and Y) should be linear cell-array(s) of matrices for which the number of rows 
% should be the same for all cells

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellrows');
end

y = cellfun(@colr, x, repmat(mat2cell(cols(:),length(cols),1),nx), 'UniformOutput', 0);

function [y] = colr(x, cols)

if ischar(cols)
  if strcmp(cols(1),'r')
    % remove n columns on the right
    cols = str2double(cols(2:end));
    cols = size(x,2)-((cols-1):-1:0);
  elseif strcmp(cols(1),'l')
    % remove n columns on the left
    cols = str2double(cols(2:end));
    cols = 1:cols;
  else
    % numeric
  end
end

x(:,cols) = [];
y         = x;