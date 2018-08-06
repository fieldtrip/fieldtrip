 function [c] = cellnancov(x, varargin)

% [C] = CELLNANCOV(X, NORMALIZEFLAG, DIM) computes the covariance, across all cells in x along 
%        the dimension dim, accounting for NaNs
% [C] = CELLNANCOV(X, Y, NORMALIZEFLAG, DIM) computes the covariance between all cells in x and y
%
% X (and Y) should be linear cell-array(s) of matrices for which the size in at 
% least one of the dimensions should be the same for all cells 

if numel(varargin)==0
  % with a single input, default parameters
  normalizeflag = 0;
  dim           = [];
  flag          = 0;
end

if numel(varargin)>=1 && iscell(varargin{1})
  % two data arguments
  y        = varargin{1};
  varargin = varargin(2:end);
else 
  y = []; 
end

if numel(varargin)>=1 && isempty(varargin{1})
  normalizeflag = 0;
elseif numel(varargin)>=1
  normalizeflag = varargin{1};
  varargin      = varargin(2:end);
end

if numel(varargin)>=1
  dim      = varargin{1};
  varargin = varargin(2:end);
else
  dim      = [];
end

if numel(varargin>=1)
  flag = varargin{1};
else
  flag = 1;
end

[scx1, scx2] = size2(x, [], 'cell');
if isempty(dim)
  if     all(scx1==scx1(1)), dim = 2;
  elseif all(scx2==scx2(1)), dim = 1; %let second dimension prevail
  else   error('no dimension to compute covariance for');
  end
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for nancov');
end

nx   = max(nx);
if isempty(y),
  n = cov(isfinite(x),[],[],0);
  if dim==1
    n = n*sum(scx1);
  elseif dim==2
    n = n*sum(scx2);
  end
  
  for k = 1:nx
    [tmp1, tmp2] = nancovc(x{k}, dim);
    if k==1
      C = tmp1;
      M = tmp2;
    else
      C = C + tmp1;
      M = M + tmp2;
    end
  end
  Mx = M;
  My = M';
else
  n = cov(isfinite(x), isfinite(y),[],[],0);
  if dim==1
    n = n*sum(scx1);
  elseif dim==2
    n = n*sum(scx2);
  end
  for k = 1:nx
    [tmp1, tmp2, tmp3] = nancovc(x{k}, y{k}, dim);
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

if normalizeflag % || n==1
  % normalize by n
  n1 = n;
  n2 = n;
else
  % normalize by n-1
  n1 = n;
  n2 = n-1;
end

if flag
  c = (C-(Mx.*My)./n1)./n2;
else
  c = C./n;
end

function [c,mx,my] = nancovc(x, y, dim)

if nargin==2
  % only a single data matrix is provided
  dim = y;
  if islogical(x), x = double(x); end

  notfinitex    = ~isfinite(x);
  x(notfinitex) = 0;

  if dim==1
    c  = x'*x;
    for k = 1:size(x,2)
    end
  elseif dim==2
    c  = x*x';
    nx = size(x,1);
    mx = repmat(sum(x,2),[1 nx]);
    for k = 1:nx
      mx(:,k) = mx(:,k) - sum(x(:,notfinitex(k,:)),2);
    end
  end
  
else
    
  if islogical(x), x = double(x); end
  if islogical(y), y = double(y); end

  notfinitex = ~isfinite(x);
  notfinitey = ~isfinite(y);

  x(notfinitex) = 0;
  y(notfinitey) = 0;

  if dim==1
    c  = x'*y;
    for k = 1:size(x,2)
    end
  elseif dim==2
    c = x*y';
    nx = size(x,1);
    ny = size(y,1);
    mx = repmat(sum(x,2),  [1 ny]);
    my = repmat(sum(y,2)', [nx 1]);
    for k = 1:ny
      mx(:,k) = mx(:,k) - sum(x(:,notfinitey(k,:)),2);
      %tmpx(:,notfinitey(k,:)) = 0;
      %mx(:,k) = sum(tmpx,2);
      %tmpx(:,notfinitey(k,:)) = x(:,notfinitey(k,:));
      %mx(:,k) = sum(x.*double(~notfinitey(k*ones(nx,1),:)&~notfinitex),2);
    end
    for k = 1:nx
      my(k,:) = my(k,:) - sum(y(:,notfinitex(k,:)),2)';
      %tmpy(:,notfinitex(k,:)) = 0;
      %my(k,:) = sum(tmpy,2);
      %tmpy(:,notfinitex(k,:)) = y(:,notfinitex(k,:));
      %my(k,:) = sum(y.*double(~notfinitey&~notfinitex(k*ones(ny,1),:)),2)';
    end
  end
end
