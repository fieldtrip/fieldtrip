function y = logn_lpdf(x,mu,sigma)
%LOGNPDF Lognormal log-probability density function (lpdf).
%
%   Y = LOGN_LPDF(X,MU,SIGMA) Returns the log of the lognormal pdf at
%   the values in X. The mean and standard deviation of log(Y) are
%   MU and SIGMA.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Default values for MU and SIGMA are 0 and 1 respectively.

% Copyright (c) 2000 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

y = -0.5 * ((log(x) - mu)./sigma).^2 - log(x .* sqrt(2*pi) .* sigma);
