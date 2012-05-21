function [w,rho] = bsscca(X, delay)

% BSSCCA computes the mixing matrix based on the canonical correlation between a signal and its lagged-one copy. It implements the algorithm described in [1]
%
% DeClercq et al 2006, IEEE Biomed Eng 2583.

if nargin<2,
  delay = 1;
end

% hmmmm we need to observe the epochs' boundaries to not create rubbish
% support cell array input

if isa(X, 'cell')
  delay = setdiff(delay(:)', 0);
  
  n  = size(X{1},1);
  mX = cellmean(X,2);
  X  = cellvecadd(X, -mX);
  
  C = cellcovshift(X,[0 delay],2,0);
  m = size(C,1)-n;
  
  A = C;
  B = C;
  
  A(1:n,1:n) = 0;
  A((n+1):end,(n+1):end) = 0;
  
  B(1:n,(n+1):end) = 0;
  B((n+1):end,1:n) = 0;
  
else
  % input is a single data matrix assumed to be a continuous stretch 
  [n,m] = size(X);
  
  % get the means
  %m   = ones(1,m-1);
  %mX  = mean(X(:,2:end),2);   % lag zero
  %mX2 = mean(X(:,1:end-1),2); % lag one
  
  % use Borga's (2001) formulation from 'a unified approach to PCA, PLS, MLR
  % and CCA'
  A = zeros(2*n);
  B = zeros(2*n);
  
  if numel(delay)==1
    Xlag = X(:,1:(m-delay));
  else
    
  end
  
  XY = X(:,delay+(1:m))*Xlag';
  XX = X(:,delay+(1:m))*X(:,delay+(1:m))';
  YY = Xlag*Xlag';
  
  %XY = (X(:,2:end)-mX*m)*(X(:,1:end-1)-mX2*m)';
  
  A(1:n,(n+1):end) = XY;
  A((n+1):end,1:n) = XY';
  B(1:n,1:n)       = XX;
  B((n+1):end,(n+1):end) = YY;
  %B(1:n,1:n)       = (X(:,2:end)-mX*m)*(X(:,2:end)-mX*m)';
  %B((n+1):end,(n+1):end) = (X(:,1:end-1)-mX2*m)*(X(:,1:end-1)-mX2*m)';
end

[w,rho]  = eig(B\A);
[~,ix]   = sort(diag(abs(rho).^2),'descend');
w        = w(1:n,ix(2:2:end))';
w        = w(1:n,1:n);
rho      = abs(rho(ix(2:2:2*n),ix(2:2:2*n))).^2;

% normalise to unit norm
% for k= 1:size(w,1)
%   w(k,:) = w(k,:)./norm(w(k,:));
% end

function [m] = cellmean(x, dim)

% [M] = CELLMEAN(X, DIM) computes the mean, across all cells in x along 
% the dimension dim.
% 
% X should be an linear cell-array of matrices for which the size in at 
% least one of the dimensions should be the same for all cells 

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellmean');
end

if nargin==1,
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; %let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   error('no dimension to compute mean for');
  end
end

nx   = max(nx);
nsmp = cellfun('size', x, dim);
ssmp = cellfun(@sum,   x, repmat({dim},1,nx), 'UniformOutput', 0);
m    = sum(cell2mat(ssmp), dim)./sum(nsmp);  


function [y] = cellvecadd(x, v)

% [Y]= CELLVECADD(X, V) - add vector to all rows or columns of each matrix 
% in cell-array X

% check once and for all to save time
persistent bsxfun_exists;
if isempty(bsxfun_exists); 
    bsxfun_exists=(exist('bsxfun')==5); 
    if ~bsxfun_exists; 
        error('bsxfun not found.');
    end
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellmean');
end

if ~iscell(v),
  v = repmat({v}, nx);
end

sx1 = cellfun('size', x, 1);
sx2 = cellfun('size', x, 2);
sv1 = cellfun('size', v, 1);
sv2 = cellfun('size', v, 2);
if all(sx1==sv1) && all(sv2==1),    
  dim = mat2cell([ones(length(sx2),1) sx2(:)]', repmat(2,nx(1),1), repmat(1,nx(2),1)); 
elseif all(sx2==sv2) && all(sv1==1),
  dim = mat2cell([sx1(:) ones(length(sx1),1)]', repmat(2,nx(1),1), repmat(1,nx(2),1));
elseif all(sv1==1) && all(sv2==1),
  dim = mat2cell([sx1(:) sx2(:)]'', nx(1), nx(2));
else   error('inconsistent input');
end  

y  = cellfun(@bsxfun, repmat({@plus}, nx), x, v, 'UniformOutput', 0);
%y = cellfun(@vplus, x, v, dim, 'UniformOutput', 0);

function [c] = cellcov(x, y, dim, flag)

% [C] = CELLCOV(X, DIM) computes the covariance, across all cells in x along 
% the dimension dim. When there are three inputs, covariance is computed between
% all cells in x and y
% 
% X (and Y) should be linear cell-array(s) of matrices for which the size in at 
% least one of the dimensions should be the same for all cells 


if nargin<4 && iscell(y)
  flag = 1;
elseif nargin<4 && isnumeric(y)
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
  error('incorrect input for cellmean');
end

if flag,
  mx   = cellmean(x, 2);
  x    = cellvecadd(x, -mx);
  if ~isempty(y),
    my = cellmean(y, 2);
    y  = cellvecadd(y, -my);
  end
end

nx   = max(nx);
nsmp = cellfun('size', x, dim);
if isempty(y), 
  csmp = cellfun(@covc, x, repmat({dim},1,nx), 'UniformOutput', 0);
else
  csmp = cellfun(@covc, x, y, repmat({dim},1,nx), 'UniformOutput', 0);
end
nc   = size(csmp{1});
c    = sum(reshape(cell2mat(csmp), [nc(1) nc(2) nx]), 3)./sum(nsmp); 

function [c] = covc(x, y, dim)

if nargin==2,
  dim = y;
  y   = x;
end

if dim==1,
  c = x'*y;
elseif dim==2,
  c = x*y';
end

function [c] = cellcovshift(x, shift, dim, flag)

% [C] = CELLCOVSHIFT(X, SHIFT, DIM) computes the covariance, across all cells
% in x along the dimension dim. 
% 
% X should be linear cell-array(s) of matrices for which the size in at 
% least one of the dimensions is be the same for all cells 

if nargin<4,
  flag = 1;
end

if nargin<3,
  dim = find(size(x{1})>1, 1, 'first');
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1) || ndims(x{1})>2,
  error('incorrect input for cellcovshift');
end

% shift the time axis
y = cellshift(x, shift, dim);

% compute covariance
c = cellcov(y, dim, flag);

function [x] = cellshift(x, shift, dim, maxshift)

% CELLSHIFT(X, SHIFT, DIM) 

if numel(shift)>1
  x = x(:)';     
  y = cell(numel(shift),numel(x));
  for k = 1:numel(shift)
    y(k,:) = cellshift(x, shift(k), dim, max(abs(shift)));
  end
  for k = 1:size(y,2)
    y{1,k} = cell2mat(y(:,k));
    for m = 2:size(y,1)
      y{m,k} = nan;
    end
  end
  x = y(1,:);
  return;
end

if nargin<4
  maxshift = max(abs(shift));
end

if any(abs(shift))>abs(maxshift)
  error('the value for maxshift should be >= shift');
end

if nargin<3
  dim = find(size(x{1})>1, 1, 'first');
end

maxshift = abs(maxshift);
if numel(maxshift)==1, maxshift = maxshift([1 1]); end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellshift');
end

n    = numel(x);
nsmp = cellfun('size', x, dim);
beg1 = ones(1,n) + shift + maxshift(1);
end1 = nsmp      + shift - maxshift(2);

switch dim
case 1
  for k = 1:n
    x{k} = x{k}(beg1(k):end1(k),:);
  end
case 2
  for k = 1:n
    x{k} = x{k}(:,beg1(k):end1(k));
  end
otherwise
  error('dimensionality of >2 is not supported');
end

