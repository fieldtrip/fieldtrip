% Copyright (C) 2012 Rik Wehbring
% Copyright (C) 1995-2012 Kurt Hornik
%
% This file is part of Octave.
%
% Octave is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% Octave is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn {Function File} {} betainv (@var{x}, @var{a}, @var{b})
% For each element of @var{x}, compute the quantile (the inverse of
% the CDF) at @var{x} of the Beta distribution with parameters @var{a}
% and @var{b}.
% @end deftypefn

% Author: KH <Kurt.Hornik@wu-wien.ac.at>
% Description: Quantile function of the Beta distribution

function inv = betainv (x, a, b)

if (nargin ~= 3)
  print_usage ();
end

if (~isscalar (a) || ~isscalar (b))
  [retval, x, a, b] = common_size (x, a, b);
  if (retval > 0)
    error ('betainv: X, A, and B must be of common size or scalars');
  end
end

if (iscomplex (x) || iscomplex (a) || iscomplex (b))
  error ('betainv: X, A, and B must not be complex');
end

if (isa (x, 'single') || isa (a, 'single') || isa (b, 'single'))
  inv = zeros (size (x), 'single');
else
  inv = zeros (size (x));
end

k = (x < 0) | (x > 1) | ~(a > 0) | ~(b > 0) | isnan (x);
inv(k) = NaN;

k = (x == 1) & (a > 0) & (b > 0);
inv(k) = 1;

k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0));
if (any (k))
  if (~isscalar (a) || ~isscalar (b))
    a = a(k);
    b = b(k);
    y = a ./ (a + b);
  else
    y = a / (a + b) * ones (size (k));
  end
  x = x(k);
  
  if (isa (y, 'single'))
    myeps = eps ('single');
  else
    myeps = eps;
  end
  
  l = find (y < myeps);
  if (any (l))
    y(l) = sqrt (myeps) * ones (length (l), 1);
  end
  l = find (y > 1 - myeps);
  if (any (l))
    y(l) = 1 - sqrt (myeps) * ones (length (l), 1);
  end
  
  y_old = y;
  for i = 1 : 10000
    h     = (betacdf (y_old, a, b) - x) ./ betapdf (y_old, a, b);
    y_new = y_old - h;
    ind   = find (y_new <= myeps);
    if (any (ind))
      y_new (ind) = y_old (ind) / 10;
    end
    ind = find (y_new >= 1 - myeps);
    if (any (ind))
      y_new (ind) = 1 - (1 - y_old (ind)) / 10;
    end
    h = y_old - y_new;
    if (max (abs (h)) < sqrt (myeps))
      break;
    end
    y_old = y_new;
  end
  
  inv(k) = y_new;
end

end