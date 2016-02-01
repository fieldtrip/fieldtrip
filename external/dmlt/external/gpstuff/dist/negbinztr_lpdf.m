function y = negbinztr_lpdf(x,l,r)
%NEGBINZTR_LPDF Zero trunc. negative binomial log probability density function
%
%  Description
%    Y = NEGBINZTR_LPDF(X,L,R) returns the log of zero truncated
%    Negative binomial probability density function with location
%    parameter L and dispersion parameter R (0<R<infty).
%
%    Negative binomial has different parameterizations and we use the form
%      p(x|l,r) = (r/(r+l))^r * gamma(r+y)/(gamma(r)*gamma(y+1))
%                             * (l/(r+l))^y
%    which approaches Poisson distribution when R goes to infinity.
%    Zero truncated Negative Binomial has p(x==0|l,r)=0.
%
%    The size of Y is the common size of X, L and R. A scalar input
%    functions as a constant matrix of the same size as the other
%    input.
%
%     Note that the density function is zero unless X is an
%     integer.
%
%  See also 
%    NEGBINZTR_PDF, NEGBIN_LPDF, NEGBIN_PDF

% Copyright (c) 2010 Jarno Vanhatalo
% Copyright (c) 2010-2011  Aki Vehtari

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
    y(x==1) = 0;
  elseif l>0
    k = (x >= 1 & x == round(x));
    if (any(k))
      lp0=r.*(log(r) - log(r+l));
      y(k)=lp0 + gammaln(r+x(k)) - gammaln(r) - gammaln(x(k)+1) + x(k).*(log(l) - log(r+l)) -log(1-exp(lp0));
    end
  end
else
  y = repmat(-Inf,size(x));
  y(l < 0) = NaN;
  y(x==1 & l==0) = 0;

  k = (x >= 1 & x == round(x) & l > 0);
  if (any(k))
    lp0=r(k).*(log(r(k)) - log(r(k)+l(k)));
    y(k)=lp0 + gammaln(r(k)+x(k)) - gammaln(r(k)) - gammaln(x(k)+1) + x(k).*(log(l(k)) - log(r(k)+l(k))) -log(1-exp(lp0));
  end
end

