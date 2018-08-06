function [c] = nancov(x, varargin)

% [C] = NANCOV(X, NORMALIZEFLAG, DIM) computes the covariance, across all cells in x along 
%        the dimension dim, accounting for NaNs
% [C] = NANCOV(X, Y, NORMALIZEFLAG, DIM) computes the covariance between all cells in x and y
% 
% X (and Y) should be linear cell-array(s) of matrices for which the size in at 
% least one of the dimensions should be the same for all cells. 

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
if ~iscell(x) || length(nx)>2 || all(nx>1)
  error('incorrect input for nancov');
end

nx   = max(nx);
if isempty(y)
  n = cov(isfinite(x),1,[],0);
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
  n = cov(isfinite(x), isfinite(y),1,[],0);
  if dim==1
    n = n*sum(scx1);
  elseif dim==2
    n = n*sum(scx2);
    tmp1 = zeros(size(x{1},1),size(y{1},1));
    tmp2 = tmp1;
    tmp3 = tmp1;
  end
  for k = 1:nx
    if all(~isfinite(x{k}(:)))||all(~isfinite(y{k}(:)))
      tmp1(:) = 0; tmp2(:) = 0; tmp3(:) = 0;
    else
      [tmp1, tmp2, tmp3] = nancovc(x{k}, y{k}, dim);
    end
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
  c = C./n2;
end

function [c,mx,my] = nancovc(x, y, dim)

if nargin==2
  dim = y;
  y   = x;
end

if islogical(x), x = double(x); end
if islogical(y), y = double(y); end

finitex = isfinite(x);
finitey = isfinite(y);

x(~finitex)=0;
y(~finitey)=0;

if dim==1
  error('not yet supported');
  c  = x'*y;
  for k = 1:size(x,2)
    z=1;
    
  end
elseif dim==2
  c = x*y';
  nx = size(x,1);
  ny = size(y,1);
  mx = zeros(nx,ny);
  my = zeros(nx,ny);
  
  % check whether nans are in the columns
  %nancolumnx = all(ismember(sum(finitex), [0 nx]));
  %nancolumny = all(ismember(sum(finitey), [0 ny])); 
  sx = sum(finitex);
  sy = sum(finitey);
  nancolumnx = all(sx==0|sx==nx);
  nancolumny = all(sy==0|sy==ny);
  
  if nancolumny
    mx = repmat(sum(x(:,finitey(1,:)),2),[1 ny]);
  else
    for k = 1:ny
      mx(:,k) = sum(x(:,finitey(k,:)),2);
    end
  end
  if nancolumnx
    my = repmat(sum(y(:,finitex(1,:)),2)',[nx 1]);
  else
    for k = 1:nx
      my(k,:) = sum(y(:,finitex(k,:)),2)';
      %my(k,:) = sum(y.*double(finitey&finitex(k*ones(ny,1),:)),2)';
    end
  end
end
