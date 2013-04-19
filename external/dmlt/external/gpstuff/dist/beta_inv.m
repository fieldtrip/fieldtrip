function binv = beta_inv (x, a, b)
%  BETA_INV Inverse of the beta cumulative distribution function (cdf).
%     X = BETA_INV(P,A,B) returns the inverse of the beta cdf with
%     parameters A and B at the values in P.

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
binv = zeros(n,nin);

k = find ((x<0) || (x>1) || (a<0) || (b<0) || isnan(x));
if (any(k))
  binv(k) = NaN;
end

k = find ((x==1) && (a>0) && (b>0));
if (any (k))
  binv(k) = 1;
end

k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0));
if (any (k))
  if (~isscalar(a) || ~isscalar(b))
    a = a(k);
    b = b(k);
    y = a./(a+b);
  else
    y = a/(a + b)*ones(size(k));
  end
  x = x(k);
  l = find (y < eps);
  if (any (l))
    y(l) = sqrt (eps) * ones (length (l), 1);
  end
  l = find (y > 1 - eps);
  if (any (l))
    y(l) = 1 - sqrt (eps) * ones (length (l), 1);
  end
  
  y_old = y;
  for i = 1 : 10000
    h     = (beta_cdf (y_old, a, b) - x) ./ beta_pdf (y_old, a, b);
    y_new = y_old - h;
    ind   = find (y_new <= eps);
    if (any (ind))
      y_new (ind) = y_old (ind) / 10;
    end
    ind = find (y_new >= 1 - eps);
    if (any (ind))
      y_new (ind) = 1 - (1 - y_old (ind)) / 10;
    end
    h = y_old - y_new;
    if (max (abs (h)) < sqrt (eps))
      break;
    end
    y_old = y_new;
  end
  
  binv(k) = y_new;
end

end