function tinv = t_inv (x, v, mu, sigma)
% T_INV   Inverse of Student's T cumulative distribution function (cdf).
%    X=TINV(P,V,MU,SIGMA) returns the inverse of Student's T cdf with V degrees
%    of freedom, at the values in P, with mean MU and standard deviation
%    SIGMA.
%
%    The size of X is the common size of P and V. A scalar input
%    functions as a constant matrix of the same size as the other input.

% Copyright (c) 1995-1997, 2005-2007 Kurt Hornik

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
if nargin < 4
  sigma = 1;
end
if nargin < 3
  mu = 0;
end
if (~isscalar(mu) || ~isscalar(sigma))
  if ~((size(x,1) == size(mu,1)) && (size(mu,1) == size(sigma,1)))
    error ('norm_inv: x, mu and sigma must be of common size or scalars');
  end
end


x = (x-mu)./sigma;

if (nargin < 2)
  error('must give atleast 2 input arguments')
end

if (~isscalar(v))
  
  if (size(x,1) ~= size(v,1))
    error ('tinv: x and v must be of common size or scalar');
  end
end

tinv = zeros(size(x));

k = find( (x<0) || (x > 1) || isnan(x) || (v<0));
if (any(k))
  tinv(k) = NaN;
end

k = find ((x == 0) & (v > 0));
if (any (k))
  tinv(k) = -Inf;
end

k = find ((x == 1) & (v > 0));
if (any (k))
  tinv(k) = Inf;
end

k = find ((x > 0) & (x < 1) & (v > 0) & (v < 10000));
if (any (k))
  if (isscalar (v))
    tinv(k) = (sign (x(k) - 1/2) ...
    .* sqrt (v .* (1 ./ beta_inv (2*min(x(k),1-x(k)), ...
    v/2,1/2)-1)));
  else
    tinv(k) = (sign(x(k)-1/2) ...
    .* sqrt (v(k).*(1./beta_inv(2*min(x(k), 1-x(k)), ...
    v(k)/2,1/2)-1)));
  end
end

% For large v, use the quantiles of the standard normal
k = find ((x > 0) & (x < 1) & (v >= 10000));
if (any (k))
  tinv(k) = sqrt (2)*erfinv(2*x(k)-1);
end

end

% inv = sqrt (2) * erfinv (2 * x - 1);