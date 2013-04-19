function inv = nbin_inv (x, r, p)
%NBIN_INV     Inverse of Negative binomial cumulative distribution function (inv).
%
%   Y = NBIN_CDF(X,R,P) Returns inverse of the Negative binomial cdf with
%   parameters R and P, at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

% Copyright (c) 1995-1997, 2007 Kurt Hornik

% This software is distributed under the GNU General Public
% Licence (version 3 or later); please refer to the file
% Licence.txt, included with the software, for details.

if (nargin ~= 3)
  error('3 parameters must be provided');
end

if (~isscalar(r) || ~isscalar(p))
  if ((size(x,1) ~= size(r,1)) || (size(r,1) ~= size(p,1)))
    error ('nbinpdf: x, r and p must be of common size or scalar');
  end
end

inv = zeros (size (x));

k = find(isnan(x) | (x<0) | (x>1) | (r<1) | isinf(r) | (p<0) | (p>1));
if (any(k))
  inv(k) = NaN;
end

k = find((x==1) & (r>0) & (r<Inf) & (p>=0) & (p<=1));
if (any(k))
  inv(k) = Inf;
end

k = find((x>=0) & (x<1) & (r>0) & (r<Inf) & (p>0) & (p<=1));
if (any (k))
  m = zeros (size (k));
  x = x(k);
  if (isscalar (r) && isscalar (p))
    s = p^r*ones(size(k));
    while (1)
      l = find(s<x);
      if (any(l))
        m(l) = m(l) + 1;
        s(l) = s(l) + nbin_pdf(m(l), r, p);
      else
        break;
      end
    end
  else
    r = r(k);
    p = p(k);
    s = p.^r;
    while (1)
      l = find(s < x);
      if (any(l))
        m(l) = m(l) + 1;
        s(l) = s(l) + nbin_pdf(m(l), r(l), p(l));
      else
        break;
      end
    end
  end
  inv(k) = m;
end

end