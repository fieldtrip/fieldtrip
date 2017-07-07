function [r] = mvnrnd(mu,sigma,varargin)
%MVNRND Random vectors from the multivariate normal distribution. This is an 
%   open source version that emulates a subpart of the behavior of the same
%   name function from the MATLAB stats-toolbox.
%   It emulates the 3 input, 1 output argument MATLAB-version, 
%   where the 3 input argument is the number of samples. If
%   more than three input arguments are provided, an error is thrown. Also,
%   the input argument SIGMA cannot be 3D.
%
%   R = MVNRND(MU,SIGMA) returns an N-by-D matrix R of random vectors
%   chosen from the multivariate normal distribution with mean vector MU,
%   and covariance matrix SIGMA.  MU is an N-by-D matrix, and MVNRND
%   generates each row of R using the corresponding row of MU.  SIGMA is a
%   D-by-D symmetric positive semi-definite matrix.
%   If the covariance matrix is diagonal, containing
%   variances along the diagonal and zero covariances off the diagonal,
%   SIGMA may also be specified as a 1-by-D matrix,
%   containing just the diagonal. If MU is a 1-by-D vector, MVNRND
%   replicates it to match the trailing dimension of SIGMA.

%   Copyright 2017, Jan-Mathijs Schoffelen

if nargin < 2 || isempty(mu) || isempty(sigma)
    error('at least 2 input arguments are needed');
elseif ndims(mu) > 2
    error('the first input argument has an unexpected size');
elseif ndims(sigma) > 3
    error('the second input argument has an unexpected size');
end

if numel(varargin)>1, 
    error('too many inputs');
end

nsmp  = varargin{1};
if isempty(nsmp)
    nsmp = 1;
end
[n,d] = size(mu);
siz   = [size(sigma) 1];
if siz(1)==1 && siz(2)>1
    % Just the diagonal of Sigma has been passed in.
    newsigma = zeros([siz(1) siz(1) siz(3:end)]);
    for k = 1:siz(3)
      newsigma(:,:,k) = diag(sigma(:,:,k));
    end
    sigma = newsigma; clear newsigma;
    siz   = [size(sigma) 1];
end

if siz(1)~=siz(2)
    error('the covariance matrix supplied has an unexpected size');
end

% Special case: if mu is a column vector, then use sigma to try
% to interpret it as a row vector.
if d == 1 && siz(1) == n
    mu = mu';
    [n,d] = size(mu);
end

% Single covariance matrix
if siz(3)==1
    T = chol(sigma); % origin mvnrnd uses cholcov from stats-toolbox, which is more robust
    r = randn(nsmp,size(T,1)) * T + mu(ones(nsmp,1),:);
end
    
