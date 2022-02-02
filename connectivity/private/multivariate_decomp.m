function [E, D] = multivariate_decomp(C,x,y,method,realflag,thr,fastflag)

% MULTIVARIATE_DECOMP does a linear decomposition of multivariate time series,
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
%             linear regression). In the case of mlr-like decompositions,
%             the indices for x reflect the independent variable)
%   realflag = true (default) or false. Do the operation on the real part
%              of the matrix if the input matrix is complex-valued
%   fastflag = true (default) or false. Compute the solution without an
%              eigenvalue decomposition (only when numel(x)==1)
%
% The implementation is based on Borga 2001, Canonical correlation, a
% tutorial (can be found online).
%
% Output arguments:
%   E = projection matrix (not necessarily normalized). to get the orientation,
%       do orix = E(x,1)./norm(E(x,1)), and oriy = E(y,1)./norm(E(y,1));
%   D = diagonal matrix with eigenvalues

if nargin<7
  fastflag = true;
end

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

A      = C;
A(x,x) = 0; % put the auto covariances to 0, keep the cross-block covariances
A(y,y) = 0;    

switch method
  case {'pls' 'plssvd' 'plsridge' 'plsqridge'}
    % partial least-squares
    B = eye(numel(x)+numel(y));
  case {'cca' 'ccasvd' 'ccaridge' 'ccaqridge'}
    % canonical correlation
    B = C;
    B(x,y) = 0; % put the cross-block covariances to 0, keep the auto covariances
    B(y,x) = 0;
  
  case {'mlr' 'mlrsvd' 'mlrridge' 'mlrqridge'}
    % regression
    B = zeros(size(A));
    B(x,x) = C(x,x);
    B(y,y) = eye(numel(y));
  otherwise
    error('unsupported method');
end

% regularize with a ridge
if contains(method, 'ridge')
  if contains(method, 'qridge')
    type = 'qridge';
  else
    type = 'ridge';
  end
  % regularize B with a ridge matrix
  R  = ridge_matrix(x, y, type, thr);
  B  = B+R;
end

% perform the decomposition in a svd-based subspace
if contains(method, 'svd')
  [ux,sx] = svd(B(x,x));
  [uy,sy] = svd(B(y,y));
  
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
  %U = blkdiag(ux(:,keepx),uy(:,keepy));
  U = zeros(size(B,1),sum(keepx)+sum(keepy));
  U(x,1:sum(keepx)) = ux(:,keepx);
  U(y,sum(keepx)+(1:sum(keepy))) = uy(:,keepy);
  
  A = U'*A*U;
  B = U'*B*U;
  C = U'*C*U;
  
  A = (A+A')./2;
  B = (B+B')./2;
  C = (C+C')./2;
end

if numel(y)==1 && fastflag
  % a full eigenvalue decomposition is an overkill, this solution is
  % probably correct (and much faster than eig, or even eigs(..,..,1)
  E(x,1) = (A(y,x)/B(x,x));
  E(y,1) = sqrt((A(y,x)*E(x))./B(y,y));
  D    = (B*E)\(A*E);
else
  % ad hoc check for well-behavedness of the matrix, and do the decomposition
  % the call to cond takes quite some time, here make the user responsible
  % for inputting well-behaved matrices
  %if cond(B)>1e8
  %  [E,D] = eig(pinv(B)*A);
  %else
  
  if realflag
    [E,D] = eig(real(A),real(B));%eig(real(B\A));
    E = real(E);
    D = real(D);
  else
    [E,D] = eig(A,B);%eig(B\A); % not using the backslash is faster
  end
  %end
end

n          = min(numel(x),numel(y));
n          = min(n, numel(diag(D)));
[D, order] = sort(diag(D), 'descend');
E          = E(:, order(1:n));
D          = D(1:n);

if ~contains(method, 'svd')
  [E(x,:),norm_x] = normc(E(x,:)); % norm normalise the coefficients per block
  [E(y,:),norm_y] = normc(E(y,:));
  if contains(method, 'mlr')
    % scale the weights for the independent variable, such that they reflect
    % proper beta weights
    %E(x,:) = E(x,:).*D(:)'.*(norm_x./norm_y);
    E(x,:) = E(x,:)*diag(norm_x./norm_y)*diag(D(:));
    
    % output the ssq error in D
    D = 1-diag(E(y,:)'*C(y,x)*E(x,:))./diag(E(y,:)'*C(y,y)*E(y,:));
  end
else
  indx_x = 1:sum(keepx);
  indx_y = sum(keepx)+(1:sum(keepy));
  
  [E(indx_x,:), norm_x] = normc(E(indx_x,:));
  [E(indx_y,:), norm_y] = normc(E(indx_y,:));
  
  if contains(method, 'mlr')
    % scale the weights for the independent variable, such that they reflect
    % proper beta weights
    E(x,:) = E(x,:)*diag(norm_x./norm_y)*diag(D(:));
    
    % output the ssq error in D
    D = 1-diag(E(y,:)'*C(y,x)*E(x,:))./diag(E(y,:)'*C(y,y)*E(y,:));
  end
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
    R = thr;%eye(nx+ny);
  case 'qridge'
    rix = eye(nx).*2 + diag(ones(nx-1,1).*-1,1) + diag(ones(nx-1,1).*-1,-1);
    rix(1)   = 1;
    rix(end) = 1;
    riy = eye(ny).*2 + diag(ones(ny-1,1).*-1,1) + diag(ones(ny-1,1).*-1,-1);
    riy(1)   = 1;
    riy(end) = 1;
    R        = zeros(nx+ny);
    R(x,x)   = rix;
    R(y,y)   = riy;
    R        = R*thr;
end

function [y, x_norm] = normc(x)

% NORMC column-wise norm normalisation

x_norm = sqrt(sum(x.^2));
y      = x*diag(1./x_norm);
