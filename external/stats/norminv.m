% Copyright (C) 2012 Rik Wehbring
% Copyright (C) 1995-2016 Kurt Hornik
% Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
%
% This file is part of the statistics package for GNU Octave.
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, see <http://www.gnu.org/licenses/>.
%
% -*- texinfo -*-
% @deftypefn  {statistics} {@var{x} =} norminv (@var{p})
% @deftypefnx {statistics} {@var{x} =} norminv (@var{p}, @var{mu})
% @deftypefnx {statistics} {@var{x} =} norminv (@var{p}, @var{mu}, @var{sigma})
%
% Inverse of the normal cumulative distribution function (iCDF).
%
% For each element of @var{p}, compute the quantile (the inverse of the CDF) of
% the normal distribution with mean @var{mu} and standard deviation @var{sigma}.
% The size of @var{p} is the common size of @var{p}, @var{mu} and @var{sigma}.
% A scalar input functions as a constant matrix of the same size as the other
% inputs.
%
% Default values are @var{mu} = 0, @var{sigma} = 1.
%
% The default values correspond to the standard normal distribution and
% computing its quantile function is also possible with the @code{probit}
% function, which is faster but it does not perform any input validation.
%
% Further information about the normal distribution can be found at
% @url{https://en.wikipedia.org/wiki/Normal_distribution}
%
% @seealso{norminv, normpdf, normrnd, normfit, normlike, normstat, probit}
% @end deftypefn

function x = norminv (p, mu, sigma)

  % Check for valid number of input arguments
  if nargin < 1
    error('norminv: function called with too few input arguments.');
  end

  % Add defaults (if missing input arguments)
  if nargin < 2
    mu = 0;
  end
  if nargin < 3
    sigma = 1;
  end

  % Check for common size of P, MU, and SIGMA
  if ~isscalar(p) || ~isscalar(mu) || ~isscalar (sigma)
    [retval, p, mu, sigma] = common_size(p, mu, sigma);
    if (retval > 0)
      error('norminv: P, MU, and SIGMA must be of common size or scalars.');
    end
  end

  % Check for P, MU, and SIGMA being reals
  if iscomplex(p) || iscomplex(mu) || iscomplex(sigma)
    error ('norminv: P, MU, and SIGMA must not be complex.');
  end

  % Check for class type
  if isa(p, 'single') || isa(mu, 'single') || isa(sigma, 'single')
    x = NaN(size(p), 'single');
  else
    x = NaN(size(p));
  end

  % Compute normal iCDF
  if isscalar(mu) && isscalar(sigma)
    if isfinite(mu) && sigma > 0 && sigma < Inf
      x = mu + sigma * (-sqrt (2) * erfcinv(2 * p));
    end
  else
    k = isfinite(mu) & sigma > 0 & sigma < Inf;
    x(k) = mu(k) + sigma(k) .* (-sqrt (2) * erfcinv(2 * p(k)));
  end

end

%!demo
%! ## Plot various iCDFs from the normal distribution
%! p = 0.001:0.001:0.999;
%! x1 = norminv (p, 0, 0.5);
%! x2 = norminv (p, 0, 1);
%! x3 = norminv (p, 0, 2);
%! x4 = norminv (p, -2, 0.8);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c")
%! grid on
%! ylim ([-5, 5])
%! legend ({"μ = 0, σ = 0.5", "μ = 0, σ = 1", ...
%!          "μ = 0, σ = 2", "μ = -2, σ = 0.8"}, "location", "northwest")
%! title ("Normal iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

% Test output
%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (norminv (p, ones (1,5), ones (1,5)), [NaN -Inf 1 Inf NaN])
%!assert (norminv (p, 1, ones (1,5)), [NaN -Inf 1 Inf NaN])
%!assert (norminv (p, ones (1,5), 1), [NaN -Inf 1 Inf NaN])
%!assert (norminv (p, [1 -Inf NaN Inf 1], 1), [NaN NaN NaN NaN NaN])
%!assert (norminv (p, 1, [1 0 NaN Inf 1]), [NaN NaN NaN NaN NaN])
%!assert (norminv ([p(1:2) NaN p(4:5)], 1, 1), [NaN -Inf NaN Inf NaN])
%!assert (norminv (p), probit (p))
%!assert (norminv (0.31254), probit (0.31254))

% Test class of input preserved
%!assert (norminv ([p, NaN], 1, 1), [NaN -Inf 1 Inf NaN NaN])
%!assert (norminv (single ([p, NaN]), 1, 1), single ([NaN -Inf 1 Inf NaN NaN]))
%!assert (norminv ([p, NaN], single (1), 1), single ([NaN -Inf 1 Inf NaN NaN]))
%!assert (norminv ([p, NaN], 1, single (1)), single ([NaN -Inf 1 Inf NaN NaN]))

% Test input validation
%!error<norminv: function called with too few input arguments.> norminv ()
%!error<norminv: P, MU, and SIGMA must be of common size or scalars.> ...
%! norminv (ones (3), ones (2), ones (2))
%!error<norminv: P, MU, and SIGMA must be of common size or scalars.> ...
%! norminv (ones (2), ones (3), ones (2))
%!error<norminv: P, MU, and SIGMA must be of common size or scalars.> ...
%! norminv (ones (2), ones (2), ones (3))
%!error<norminv: P, MU, and SIGMA must not be complex.> norminv (i, 2, 2)
%!error<norminv: P, MU, and SIGMA must not be complex.> norminv (2, i, 2)
%!error<norminv: P, MU, and SIGMA must not be complex.> norminv (2, 2, i)
