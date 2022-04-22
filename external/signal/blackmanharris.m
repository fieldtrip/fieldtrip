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
% @deftypefn {Function File} {[@var{w}] =} blackmanharris(@var{L})
% Compute the Blackman-Harris window.
% @seealso{rectwin,  bartlett}
% @end deftypefn

function [w] = blackmanharris(L, opt)
  if (nargin < 1)
    help(mfilename);
  elseif(~ isscalar(L))
    error('L must be a number');
  end % if
  
  N = L - 1;
  if (nargin == 2)
    switch (opt)
      case 'periodic'
        N = L;
      case 'symmetric'
        N = L - 1;
      otherwise
        error('blackmanharris: window type must be either ''periodic'' or ''symmetric''');
    end %switch
  end %if

  a0 = 0.35875;
  a1 = 0.48829;
  a2 = 0.14128;
  a3 = 0.01168;
  n = (0:(L-1))';
  
  if N==0
    w = 1;
  else
    w = a0 - a1.*cos(2.*pi.*n./N) + a2.*cos(4.*pi.*n./N) - a3.*cos(6.*pi.*n./N);
  end
end
