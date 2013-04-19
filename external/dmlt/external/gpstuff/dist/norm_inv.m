function ninv = norm_inv(x,mu,sigma)
%NORM_INV Inverse of the normal cumulative distribution function (cdf).
%   Y = NORM_INV(X,MU,SIGMA) Returns inverse of the normal cdf with
%   mean, MU, and standard deviation, SIGMA, at the values in X.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1 respectively.

% Copyright (c) 1995-1997, 2005-2007 Kurt Hornik

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.


if nargin < 3
  sigma = 1;
end
if nargin < 2
  mu = 0;
end

if (~isscalar(mu) || ~isscalar(sigma))
  if ~((size(x,1) == size(mu,1)) && (size(mu,1) == size(sigma,1)))
    error ('norm_inv: x, mu and sigma must be of common size or scalars');
  end
end

[n,nin] = size (x);
ninv = zeros (n,nin);

if (isscalar (mu) && isscalar (sigma))
  if (find (isinf (mu) || isnan (mu) || (sigma < 0) || isinf(sigma)))
    ninv = NaN*ones(n,nin);
  else
    ninv =  mu+sigma.*sqrt(2)*erfinv(2*x-1);
  end
else
  k = find(isinf(mu) || isnan(mu) || (sigma < 0) || isinf(sigma));
  if (any(k))
    ninv(k) = NaN;
  end
  
  k = find(~isinf(mu) && ~isnan(mu) && (sigma > 0) && (sigma < Inf));
  if (any (k))
    ninv(k) = mu(k)+sigma(k).*sqrt(2)*erfinv(2*x(k)-1);
  end
end

k = find ((sigma == 0) && (x > 0) && (x < 1));
if (any(k))
  ninv(k) = mu(k);
end

ninv((sigma == 0) & (x == 0)) = -Inf;
ninv((sigma == 0) & (x == 1)) = Inf;

end