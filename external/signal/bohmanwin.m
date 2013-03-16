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
% @deftypefn {Function File} {[@var{w}] =} bohmanwin(@var{L})
% Compute the Bohman window of lenght L.
% @seealso{rectwin,  bartlett}
% @end deftypefn

function [w] = bohmanwin(L)
  if (nargin < 1)
    help(mfilename)
  elseif(~ isscalar(L))
    error('L must be a number');
  elseif(L < 0)
    error('L must be positive');
  end

  if(L ~= floor(L))
    L = round(L);
    warning('L rounded to the nearest integer.');
  end

  if(L == 0)
    w = [];

  elseif(L == 1)
    w = 1;

  else
    N = L-1;
    n = -N/2:N/2;

    w = (1-2.*abs(n)./N).*cos(2.*pi.*abs(n)./N) + (1./pi).*sin(2.*pi.*abs(n)./N);
    w(1) = 0;
    w(length(w))=0;
    w = w';
  end
end
