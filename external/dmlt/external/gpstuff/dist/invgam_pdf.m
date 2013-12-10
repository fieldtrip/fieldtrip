function y = invgam_pdf(x,a,b)
%INVGAM_PDF Inverse-Gamma probability density function.
%
%   Y = INVGAM_PDF(X,A,B) returns the inverse- gamma probability
%   density function with parameters A and B, at the values in X.
%
%   Note: Parameterization as in (Neal, 1996).
%      A is mean of the distribution
%      B is degrees of freedom

% Copyright (c) 2000 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if nargin < 2, 
   error('Requires at least two input arguments.'); 
end

y = (1-b./2) .* log(x) - (a.*b./2./x) - (b/2) .* log(2./(a.*b)) -gammaln(b/2);
y=exp(y);
