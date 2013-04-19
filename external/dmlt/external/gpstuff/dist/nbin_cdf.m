function cdf = nbin_cdf (x, r, p)
%NBIN_CDF     Negative binomial cumulative distribution function (cdf).
%
%   Y = NBIN_CDF(X,R,P) Returns the Negative binomial cdf with
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
    error ('nbincdf: x, r and p must be of common size or scalar');
  end
end

cdf = zeros (size (x));

k = find (isnan (x) | (r < 1) | (r == Inf) | (p < 0) | (p > 1));
if (any (k))
  cdf(k) = NaN;
end

k = find (isinf(x) & (r > 0) & (r < Inf) & (p >= 0) & (p <= 1));
if (any (k))
  cdf(k) = 1;
end

k = find ((x >= 0) & (x < Inf) & (x==round(x)) & (r > 0) & ~isinf(r) & (p>0) & (p<=1));
if (any (k))
  m = zeros(size(k));
  x = floor (x(k));
  y = cdf(k);
  if (isscalar(r) && isscalar(p))
    while (1)
      l = find(m <= x);
      if (any (l))
        y(l) = y(l) + nbin_pdf(m(l), r, p);
        m(l) = m(l) + 1;
      else
        break;
      end
    end
  else
    r = r(k);
    p = p(k);
    while (1)
      l = find (m <= x);
      if (any (l))
        y(l) = y(l) + nbin_pdf(m(l), r(l), p(l));
        m(l) = m(l) + 1;
      else
        break;
      end
    end
  end
  cdf(k) = y;
end

end