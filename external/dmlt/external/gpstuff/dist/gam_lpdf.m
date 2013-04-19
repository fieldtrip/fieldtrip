function y = gam_lpdf(x,a,b)
%GAM_LPDF Log of Gamma probability density function (lpdf).
%
%   Y = GAM_LPDF(X,A,B) Returns the gamma lpdf with
%   shape A and inverse scale B, at the values in X.
%
%   The size of X is the common size of the input arguments. A
%   scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   The parameterization is as in Gelman et al. (2004). Bayesian Data
%   Analysis (second edition)

% Copyright (c) 1998-2004 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if nargin < 3, 
  error('Requires three input arguments.');
end

%y = b.^a/gamma(a)*x^(a-1)*exp(-b*x);
y = a.*log(b)-gammaln(a)+(a-1).*log(x)-b.*x;
