function pdf = nbin_pdf (x, r, p)
%NBIN_PDF    Negative binomial probability density function (pdf).
%
%   Y = NBIN_PDF(X,R,P) Returns the Negative binomial pdf with
%   parameters R and P, at the values in X.
%
%   The size of Y is the common size of the input arguments. A
%   scalar input functions as a constant matrix of the same size as
%   the other inputs.
%

% Copyright (c) 1995-1997,2007 Kurt Hornik

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
% p=1-p
if (nargin ~= 3)
  error('3 parameters must be provided');
end

if (~isscalar(r) || ~isscalar(p))
  if ((size(x,1) ~= size(r,1)) || (size(r,1) ~= size(p,1)))
    error ('nbinpdf: x, r and p must be of common size or scalar');
  end
end

pdf = zeros (size (x));

k = find (isnan (x) | (r < 1) | isinf(r) | (p < 0) | (p > 1));
if (any (k))
  pdf(k) = NaN;
end

% Just for the fun of it ...
k = find (isinf(x) & (r > 0) & (r < Inf) & (p == 0));
if (any (k))
  pdf(k) = 1;
end

k = find ((x >= 0) & (x < Inf) & (x == round (x)) & (r > 0) & (r < Inf) & (p > 0) & (p <= 1));
if (any (k))
  if (isscalar (r) && isscalar (p))
    pdf(k) = p^r.*gamma(r+x(k))./(gamma(r).*gamma(x(k)+1)) .* (1-p).^x(k);
%     pdf(k) = bincoeff (-n, x(k)) .* (p ^ n) .* ((p - 1) .^ x(k));
  else
    pdf(k) = p(k).^r(k).*gamma(r(k)+x(k))./(gamma(r(k)).*gamma(x(k)+1)) .* (1-p(k)).^x(k);
%     pdf(k) = bincoeff (-n(k), x(k)) .* (p(k) .^ n(k)) .* ((p(k) - 1) .^ x(k));
  end
end

end