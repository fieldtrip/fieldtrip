function [y] = cellcrossproduct(x)

% [Y] = CELLSUBSELECT(X, BOOLVEC, DIM) outputs a cell-arry Y with the same dimension as X
% but for each input cell the n(n+1)/2 cross products of the rows are computed.
% the cross-products are ordered along the lower-triangle of the pairwise
% matrix.

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellcrossproduct');
end

n = cellfun('size', x, 1);
if ~all(n==n(1))
  error('each cell should have the same number of rows');
end
%[ix{1}(:,1),ix{1}(:,2)] = find(tril(ones(n(1))));

%y = cellfun(@crossprod, x, repmat(ix, nx), 'UniformOutput', 0);
y = cellfun(@crossprod, x, 'UniformOutput', 0);

function [y] = crossprod(x)

% FIXME only works in case of dim=2
n    = size(x);
indx = tril(ones(n(1)))==1;
y = zeros(0.5*n(1)*(n(1)+1), n(2));
for k = 1:n(2)
  tmp = tril(x(:,k)*x(:,k)');
  y(:,k) = tmp(indx);
end


% function [y] = crossprod(x, ix)
% 
% %FIXME works only in case of dim=2
% n = size(x);
% y = zeros(0.5*n(1)*(n(1)+1), n(2));
% for k = 1:size(ix,1)
%   y(k,:) = x(ix(k,1),:).*x(ix(k,2),:);
% end

