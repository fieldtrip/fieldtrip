function [c] = cellcov(x, y, dim, flag)

% CELLCOV computes the covariance, across all cells in x along 
% the dimension dim. When there are three inputs, covariance is computed between
% all cells in x and y
% 
% X (and Y) should be linear cell-array(s) of matrices for which the size in at 
% least one of the dimensions should be the same for all cells 


if nargin<4 && iscell(y)
  flag = 1;
elseif nargin==3 && isnumeric(y)
  flag = dim;
end

if nargin<3 && iscell(y)
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; %let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   error('no dimension to compute covariance for');
  end
elseif nargin<=3 && isnumeric(y)
  dim = y;
end

if isnumeric(y), y = []; end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellcov');
end

nx   = max(nx);
nsmp = cellfun('size', x, dim);
n    = sum(nsmp);
if isempty(y),  
  for k = 1:nx
    [tmp1, tmp2] = covc(x{k}, dim);
    if k==1
      C = tmp1;
      M = tmp2;
    else
      C = C + tmp1;
      M = M + tmp2;
    end
  end
  Mx = M;
  My = M;
else
  for k = 1:nx
    [tmp1, tmp2, tmp3] = covc(x{k}, y{k}, dim);
    if k==1
      C  = tmp1;
      Mx = tmp2;
      My = tmp3;
    else  
      C  = C  + tmp1;
      Mx = Mx + tmp2;
      My = My + tmp3;
    end  
  end
end

if flag
  c = (C-(Mx(:)*My(:)')./n)./n;
else
  c = C./n;
end
