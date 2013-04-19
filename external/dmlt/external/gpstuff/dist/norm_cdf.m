function y = norm_cdf(x,mu,sigma)
%NORM_CDF Normal cumulative probability density function (cdf).
%
%   Y = NORMCDF(X,MU,SIGMA) Returns the normal cdf with
%   mean, MU, and standard deviation, SIGMA, at the values in X.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Default values for MU and SIGMA are 0 and 1 respectively.

% Copyright (c) 1998-2004,2011 Aki Vehtari

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

y = 0.5 * erfc(-((x-mu) ./ sigma) ./ sqrt(2));
