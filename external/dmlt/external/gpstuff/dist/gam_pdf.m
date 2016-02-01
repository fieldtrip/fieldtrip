function y = gam_pdf(x,s2,nu)
%GAM_PDF Gamma probability density function (pdf).
%
%   Y = GAM_PDF(X,S2,NU) Returns the gamma pdf with
%   X ~ Gamma(s2, nu), where s2 is the shape and nu the inverse scale.
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
y = exp(s2.*log(nu)-gammaln(s2)+(s2-1).*log(x)-nu.*x);
