function tcdf = t_cdf (x, v, mu, sigma)
%T_PDF     Student's t cumulative distribution function (cdf).
%
%   Y = T_PDF(X,V,MU,SIGMA) Returns the Student's t pdf with
%   V degrees of freedom, MU mean and SIGMA standard deviation at X.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Default values for MU and SIGMA are 0 and 1 respectively.

% Copyright (c) 1995-1997, 2005-2007 Kurt Hornik

% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.

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

if ~isscalar(v)
  if (size(v,1) ~= size(x,1))
    error('tcdf: x and v must be of common size or scalar');
  end
end

tcdf = zeros(size (x));

k = find(isnan (x) | (v < 0));
if (any (k))
  tcdf(k) = NaN;
end

k = find ((x == Inf) & (v > 0));
if (any (k))
  cdf(k) = 1;
end

k = find ((x > -Inf) & (x < Inf) & (v > 0));
if (any (k))
  if (isscalar (v))
    tcdf(k) = betainc(1./(1+x(k).^2./v),v/2,1/2)/2;
  else
    tcdf(k) = betainc(1./(1+x(k).^2./v(k)),v(k)/2,1/2)/2;
  end
  ind = find (x(k) > 0);
  if (any (ind))
    tcdf(k(ind)) = 1 - tcdf(k(ind));
  end
end

end