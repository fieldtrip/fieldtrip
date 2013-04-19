function y = negbin_lpdf(x,l,r)
%NEGBIN_LPDF Negative binomial log probability density function
%
%   Description
%   Y = NEGBIN_LPDF(X,L,R) returns the log of Negative binomial
%   probability density function with location parameter L and
%   dispersion parameter R (0<R<infty).
%
%   Negative binomial has different parameterizations and we use the form
%    p(x|l,r) = (r/(r+l))^r * gamma(r+y)/(gamma(r)*gamma(y+1))
%                             * (l/(r+l))^y
%   which approaches Poisson distribution when R goes to infinity.
%
%   The size of Y is the common size of X, L and R. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Note that the density function is zero unless X is an integer.
%
%   See also NEGBIN_PDF

% Copyright (c) 2010 Jarno Vanhatalo
% Copyright (c) 2010-2011 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if isscalar(r) && isscalar(l)
  if l<0
    y = repmat(NaN,size(x));
  else
    y = repmat(-Inf,size(x));
  end
  if l==0
    y(x==0) = 0;
  elseif l>0
    k = (x >= 0 & x == round(x));
    if (any(k))
      y(k)=r.*(log(r) - log(r+l)) + gammaln(r+x(k)) - gammaln(r) - gammaln(x(k)+1) + x(k).*(log(l) - log(r+l));
    end
  end
else
  y = repmat(-Inf,size(x));
  y(l < 0) = NaN;
  y(x==0 & l==0) = 0;

  k = (x >= 0 & x == round(x) & l > 0);
  if (any(k))
    y(k)=r(k).*(log(r(k)) - log(r(k)+l(k))) + gammaln(r(k)+x(k)) - gammaln(r(k)) - gammaln(x(k)+1) + x(k).*(log(l(k)) - log(r(k)+l(k)));
  end
end
