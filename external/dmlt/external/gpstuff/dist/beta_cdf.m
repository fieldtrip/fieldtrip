function bcdf = beta_cdf (x, a, b)
%  BETACDF Beta cumulative distribution function.
%     P = BETACDF(X,A,B) returns the beta cumulative distribution
%     function with parameters A and B at the values in X.

% Copyright (c) 1995-1997, 2005-2007 Kurt Hornik

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
if (nargin ~= 3)
  error('must provide 3 parameters')
end

if (~isscalar (a) || ~isscalar(b))
  if ((size(a,1) ~= size(b,1)) || (size(b,1) ~= size(x,1)))
    error ('betainv: x, a and b must be of common size or scalars');
  end
end

[n,nin] = size(x);
bcdf = zeros(n,nin);

k = find ((a < 0) || (b < 0) || isnan (x));
if (any(k))
  bcdf(k) = NaN;
end

k = find((x >= 1) && (a > 0) && (b > 0));
if (any(k))
  bcdf(k) = 1;
end

k = find ((x > 0) && (x < 1) && (a > 0) && (b > 0));
if (any (k))
  if (isscalar(a) && isscalar(b))
    bcdf(k) = betainc(x(k), a, b);
  else
    bcdf(k) = betainc(x(k), a(k), b(k));
  end
end

end