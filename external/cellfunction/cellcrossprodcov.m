function [c] = cellcrossprodcov(x, p)

% CELLPRODCOV(X) computes the covariance, across all cells in x along 
% the second dimension, where x has been expanded to contain all pairwise
% cross-products (i.e. between the rows).

if nargin<2
  p = 1;
end

scx1 = cellfun('size', x, 1);
% scx2 = cellfun('size', x, 2);
if ~all(scx1==scx1(1))
  error('cellcrossprodcov operates on the second dimensions of each cell. input is incompatible');
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellcov');
end

nx   = max(nx);
nsmp = cellfun('size', x, 2);
n    = sum(nsmp);
for k = 1:nx
  [tmp1, tmp2] = crossprodcovc(x{k}, p);
  if k==1
    C = tmp1;
    M = tmp2;
  else
    C = C + tmp1;
    M = M + tmp2;
  end
end
% c = (C-(M(:)*M(:)')./n)./n;
c = C./n;