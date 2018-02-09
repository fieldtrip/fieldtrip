function [E, D] = multivariate_decomp(C,x,y,method,realflag,thr)

% MULTIVARIATE_DECOMP does a decomposition of multivariate time series,
% based on the covariance matrix.
%
% Use as:
%  [E, D] = multivariate_decomp(C,x,y,method)
%
% Input arguments:
%   C = covariance matrix (or csd) between input time series 
%   x = list of indices corresponding to group 1 
%   y = list of indices corresponding to group 2 
%   method = 'cca', or 'pls', 'mlr', decomposition method
%            (canonical correlation partial least squares, or multivariate
%             linear regression).
%   realflag = true (default) or false. Do the operation on the real part
%              of the matrix if the input matrix is complex-valued
%
% The implementation is based on Borga 2001, Canonical correlation, a
% tutorial (can be found online).
%
% Output arguments:
%   E = projection matrix (unnormalized). to get the orientation,
%       do orix = E(x,1)./norm(E(x,1)), and oriy = E(y,1)./norm(E(y,1));
%   D = diagonal matrix with eigenvalues

if nargin<6
  thr = 0.9; % percentage of variance to account for when doing ccasvd
end

if numel(thr)==1
  thr = [thr thr];
end

if nargin<5
  realflag = true;
end

if nargin<4
  method = 'cca';
end

% if realflag
%   C = real(C);
% end

A = C;
A(x,x) = 0;
A(y,y) = 0;
    
switch method
  case {'pls' 'plssvd' 'plsridge' 'plsqridge'}
    B = eye(numel(x)+numel(y));
  case {'cca' 'ccasvd' 'ccaridge' 'ccaqridge'}
    B = C;
    B(x,y) = 0;
    B(y,x) = 0;
  
  case {'mlr' 'mlrsvd' 'mlrridge' 'mlrqridge'}
    B = zeros(size(A));
    B(x,x) = C(x,x);
    B(y,y) = eye(numel(y));
  otherwise
    error('unsupported method');
end

if ~isempty(strfind(method, 'ridge'))
  if ~isempty(strfind(method, 'qridge'))
    type = 'qridge';
  else
    type = 'ridge';
  end
  % regularize B with a ridge matrix
  R  = ridge_matrix(x, y, type, thr);
  B  = B+R;
end

if ~isempty(strfind(method, 'svd'))
  [ux,sx,~] = svd(B(x,x));
  [uy,sy,~] = svd(B(y,y));
  
  if thr(1)<1
    keepx = cumsum(diag(sx))./sum(diag(sx))<=thr(1); keepx(1) = true;
  else
    keepx = false(numel(x),1);
    keepx(1:min(thr(1),numel(x))) = true;
  end
  if thr(2)<1
    keepy = cumsum(diag(sy))./sum(diag(sy))<=thr(2); keepy(1) = true;
  else
    keepy = false(numel(y),1);
    keepy(1:min(thr(2),numel(y))) = true;
  end
  U = blkdiag(ux(:,keepx),uy(:,keepy));
  
  A = U'*A*U;
  B = U'*B*U;
  
  A = (A+A')./2;
  B = (B+B')./2;
end

% ad hoc check for well-behavedness of the matrix
if cond(B)>1e8
  [E,D] = eig(pinv(B)*A);
else
  if realflag
   [E,D] = eig(real(B\A));
   E = real(E);
   D = real(D);
  else
   [E,D] = eig(B\A);
  end
end

n          = min(numel(x),numel(y));
n          = min(n, numel(diag(D)));
[D, order] = sort(diag(D), 'descend');
E          = E(:, order(1:n));
D          = D(1:n);

if ~isempty(strfind(method, 'svd'))
  E = U*E;
end

function R = ridge_matrix(x, y, type, thr)

if numel(thr)==numel(x)+numel(y)
  thr = diag(thr);
elseif numel(thr)==2
  tmp = zeros(numel(x)+numel(y),1);
  tmp(x) = thr(1);
  tmp(y) = thr(2);
  thr    = diag(tmp);
end

nx = numel(x);
ny = numel(y);
switch type
  case 'ridge'
    R = eye(nx+ny);
  case 'qridge'
    rix = eye(nx).*2 + diag(ones(nx-1,1).*-1,1) + diag(ones(nx-1,1).*-1,-1);
    rix(1)   = 1;
    rix(end) = 1;
    riy = eye(ny).*2 + diag(ones(ny-1,1).*-1,1) + diag(ones(ny-1,1).*-1,-1);
    riy(1)   = 1;
    riy(end) = 1;
    R = blkdiag(rix,riy);
end
R = R*thr;

