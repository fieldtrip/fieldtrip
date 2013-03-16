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
% @deftypefn {Function File} {[@var{w}] =} nuttallwin(@var{L})
% Compute the Blackman-Harris window defined by Nuttall of length L.
% @seealso{blackman, blackmanharris}
% @end deftypefn

function [w] = nuttallwin(L)
  if (nargin ~= 1); help(mfilename); end

  if(L < 0)
    error('L must be positive');
  end

  if(L ~= floor(L))
    L = round(L);
    warning('L rounded to the nearest integer.');
  end

  N = L-1;
  a0 = 0.355768;
  a1 = 0.487396;
  a2 = 0.144232;
  a3 = 0.012604;
  n = -N/2:N/2;
  w = a0 + a1.*cos(2.*pi.*n./N) + a2.*cos(4.*pi.*n./N) + a3.*cos(6.*pi.*n./N);
  w = w';
end
