function [c] = cov(x, varargin)

% [C] = COV(X, NORMALIZEFLAG, DIM) computes the covariance, across all cells in x along 
%        the dimension dim. W
% [C] = COV(X, Y, NORMALIZEFLAG, DIM) computes the covariance between all cells in x and y
% 
% X (and Y) should be linear cell-array(s) of matrices for which the size in at 
% least one of the dimensions should be the same for all cells 

if numel(varargin)==0
  normalizeflag = 0;
  dim           = [];
  flag          = 0;
end

if numel(varargin)>=1 && iscell(varargin{1})
  y        = varargin{1};
  varargin = varargin(2:end);
else 
  y = []; 
end

if numel(varargin)>=1 && isempty(varargin{1})
  normalizeflag = 0;
  varargin      = varargin(2:end);
elseif numel(varargin)>=1
  normalizeflag = varargin{1};
  varargin      = varargin(2:end);
end

if numel(varargin)>=1
  dim      = varargin{1};
  varargin = varargin(2:end);
end

if numel(varargin>=1)
  flag = varargin{1};
else
  flag = 1;
end

if isempty(dim)
  [scx1, scx2] = size2(x, [], 'cell');
  if     all(scx1==scx1(1)), dim = 2;
  elseif all(scx2==scx2(1)), dim = 1; %let second dimension prevail
  else   error('no dimension to compute covariance for');
  end
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1)
  error('incorrect input for cellcov');
end

nx   = max(nx);
nsmp = cellfun('size', x, dim);
n    = sum(nsmp);
if isempty(y)  
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

if normalizeflag || n==1
  % normalize by n
  n1 = n;
  n2 = n;
else
  % normalize by n-1
  n1 = n;
  n2 = n-1;
end

if flag
  c = (C-(Mx(:)*My(:)')./n1)./n2;
else
  c = C./n2;
end

function [c,mx,my] = covc(x, y, dim)

if nargin==2
  dim = y;
  y   = x;
end

if islogical(x), x = double(x); end
if islogical(y), y = double(y); end

if dim==1
  c  = x'*y;
elseif dim==2
  c = x*y';
end

mx = sum(x,dim);
my = sum(y,dim);

