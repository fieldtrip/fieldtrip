function y = poiss_lpdf(x,l)
%POISS_LPDF Poisson log-probability density function.
%
%   Description
%   Y = POISS_LPDF(X,L) returns the log of Poisson probability density 
%   function with location parameter L at the values in X.
%
%   The size of Y is the common size of X and L. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Note that the density function is zero unless X is an integer.
%
%   See also POISS_PDF

% Copyright (c) 1998-2004 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

y = repmat(-Inf,size(x));
y(l < 0) = NaN;
y(x==0 & l==0) = 0;

k = (x >= 0 & x == round(x) & l > 0);
if (any(k))
  y(k) = -l(k) +x(k).*log(l(k)) -gammaln(x(k)+1);
end
