function y = laplace_lpdf(x,mu,sigma)
%LAPLACE_LPDF    Laplace log-probability density function (lpdf).
%
%   Y = LAPLACE_LPDF(X,MU,SIGMA) Returns the log of the Laplace pdf with
%   mean, MU, and scale SIGMA, at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Default values for MU and SIGMA are 0 and 1 respectively.

% Copyright (c) 2005 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if nargin < 3, 
  sigma = 1;
end

if nargin < 2;
  mu = 0;
end

if nargin < 1, 
  error('Requires at least one input argument.');
end

y= -abs(x-mu)./sigma - log(2.*sigma);
