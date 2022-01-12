% Copyright (C) 2007 Sylvain Pelissier <sylvain.pelissier@gmail.com>
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
% @deftypefn {Function File} {[@var{w}] =} barthannwin(@var{L})
% Compute the modified Bartlett-Hann window of lenght L.
% @seealso{rectwin,  bartlett}
% @end deftypefn

function [w] = barthannwin(L)
  if (nargin < 1)
    help(mfilename);
  elseif (~ isscalar(L) || L < 0)
    error('L must be a positive integer');
  end % if
  L = round(L);
  N = L-1;
  n = 0:N;

  if N==0
    w = 1;
  else
    w = 0.62 -0.48.*abs(n./(L-1) - 0.5)+0.38*cos(2.*pi*(n./(L-1)-0.5));
    w = w';
  end
end
