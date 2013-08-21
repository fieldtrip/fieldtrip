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
% @deftypefn {Function File} {} finv (@var{x}, @var{m}, @var{n})
% For each element of @var{x}, compute the quantile (the inverse of
% the CDF) at @var{x} of the F distribution with @var{m} and @var{n}
% degrees of freedom.
% @end deftypefn

% Author: KH <Kurt.Hornik@wu-wien.ac.at>
% Description: Quantile function of the F distribution

function inv = finv (x, m, n)

if (nargin ~= 3)
  print_usage ();
end

if (~isscalar (m) || ~isscalar (n))
  [retval, x, m, n] = common_size (x, m, n);
  if (retval > 0)
    error ('finv: X, M, and N must be of common size or scalars');
  end
end

if ~(isreal (x) || isreal (m) || isreal (n))
  error ('finv: X, M, and N must not be complex');
end

if (isa (x, 'single') || isa (m, 'single') || isa (n, 'single'))
  inv = NaN (size (x), 'single');
else
  inv = NaN (size (x));
end

k = (x == 1) & (m > 0) & (m < Inf) & (n > 0) & (n < Inf);
inv(k) = Inf;

k = (x >= 0) & (x < 1) & (m > 0) & (m < Inf) & (n > 0) & (n < Inf);
if (isscalar (m) && isscalar (n))
  inv(k) = ((1 ./ betainv (1 - x(k), n/2, m/2) - 1) * n / m);
else
  inv(k) = ((1 ./ betainv (1 - x(k), n(k)/2, m(k)/2) - 1)...
    .* n(k) ./ m(k));
end

end


