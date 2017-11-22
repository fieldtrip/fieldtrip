function [unmixing, mixing, rho, compdata, time] = bsscca(X, varargin)

% BSSCCA computes the unmixing matrix based on the canonical correlation between 
% two sets of (possibly multivariate) signals, the sets may contain time shifted copies. 
% In its default, it implements the algorithm described in [1], computing the
% canonical correlation between a set of signals and their lag-one-shifted
% copy. Alternatively, if the input contains a reference signal (possibly multivariate),
% the canonical correlation between the data in X and the reference signal is computed.
% It requires JM's cellfunction toolbox on the MATLAB path:
%  (github.com/schoffelen/cellfunction.git)
%
% [1] DeClercq et al 2006, IEEE Biomed Eng 2583.

% Copyright (C) 2017, Jan-Mathijs Schoffelen

if isa(X, 'cell')
  n = size(X{1},1);
else
  n = size(X,1);
end

refdelay  = ft_getopt(varargin, 'delay',   1);
refdelay  = ft_getopt(varargin, 'refdelay', refdelay); % allow for old naming 'delay'
chandelay = ft_getopt(varargin, 'chandelay', 0);
time      = ft_getopt(varargin, 'time');
Y         = ft_getopt(varargin, 'refdata', {});
dowhiten  = istrue(ft_getopt(varargin, 'prewhiten',  0));
tol       = ft_getopt(varargin, 'tol', 1e-6);

hasrefdata = ~isempty(Y);

% hmmmm we need to observe the epochs' boundaries to not create rubbish
% support cell array input

if isa(X, 'cell')
  
  % demean
  X  = cellvecadd(X, -cellmean(X, 2));
  
  maxshift = [abs(min(min(refdelay),min(chandelay))) abs(max(max(refdelay),max(chandelay)))];
  
  if ~hasrefdata
    Y  = cellshift(X,  refdelay,  2, maxshift);
  else
    Y  = cellvecadd(Y, -cellmean(Y, 2));
    Y  = cellshift(Y, refdelay, 2, maxshift);
  end
  X  = cellshift(X, chandelay, 2, maxshift);
  
  if ~isempty(time)
    time = cellshift(time, 0, 2, maxshift);
  end
  
  if dowhiten
    % do an svd based whitening based on the relative tolerance parameter tol
    % (relative to the largest singular value)
    [ux,sx,vx] = svd(cov(X,1,2,1));
    [uy,sy,vy] = svd(cov(Y,1,2,1));
    
    sx = diag(sx);
    sy = diag(sy);
    
    keepx = sx./sx(1)>tol;
    keepy = sy./sy(1)>tol;
    
    fprintf('whitening the data, keeping %d signal components in X\n', sum(keepx));
    fprintf('whitening the data, keeping %d signal components in Y\n', sum(keepy));
    
    whiten_x = diag(1./sqrt(sx(keepx)))*ux(:,keepx)';
    whiten_y = diag(1./sqrt(sy(keepy)))*uy(:,keepy)';
    
    unwhiten_x = ux(:,keepx)*diag(sqrt(sx(keepx)));
    unwhiten_y = uy(:,keepy)*diag(sqrt(sy(keepy)));
    
    X = whiten_x*X;
    Y = whiten_y*Y;
  end

  
  % compute covariance
  C  = cov(cellcat(1,X,Y), 1, 2, 1);
  ix = 1:size(X{1},1);
  iy = ix(end)+(1:size(Y{1},1));
  
  
  XY = C(ix,iy); % cross terms XY
  XX = C(ix,ix); % auto terms X;
  YY = C(iy,iy); % auto terms Y
  YX = C(iy,ix);
  
else
  ft_error('this does not work at the moment');
  ft_warning('Running bsscca with concatenated trial in the input, represented as a single matrix, is not optimal. Consider using cellmode');
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
  
  XY = X(:,delay+(1:(m-delay)))*Xlag';
  XX = X(:,delay+(1:(m-delay)))*X(:,delay+(1:(m-delay)))';
  YY = Xlag*Xlag';
  
  A(1:n,(n+1):end) = XY;
  A((n+1):end,1:n) = XY';
  B(1:n,1:n)       = XX;
  B((n+1):end,(n+1):end) = YY;
end

XXXY = XX\XY;
YYYX = YY\YX;

[wx,rho]  = eig(XXXY*YYYX);
[wy,rho2] = eig(YYYX*XXXY);

rho      = sqrt(real(diag(rho)));
rho2     = sqrt(real(diag(rho2)));
[rho,srt1]  = sort(rho, 'descend');
[rho2,srt2] = sort(rho2, 'descend');

n  = min(numel(srt1),numel(srt2));
wx = wx(:,srt1(1:n));
wy = wy(:,srt2(1:n));
rho = rho(1:n);

% get the tranpose
wx = wx';
wy = wy';

% unmix the data
x = wx*X;
y = wy*Y;


ax = (XX*wx')/cov(x,1,2,1); % by construction wx*ax = I
ay = (YY*wy')/cov(y,1,2,1);


if dowhiten
  wx = wx*whiten_x;
  wy = wy*whiten_y;
  ax = unwhiten_x*ax;
  ay = unwhiten_y*ay;
end

if hasrefdata
  unmixing = blkdiag(wx,wy);
  mixing   = blkdiag(ax,ay);
  compdata = cellcat(1, x, y);
else
  unmixing = wx;
  mixing   = ax;
  compdata = x;
end 
