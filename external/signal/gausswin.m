% Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
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
% usage: w = gausswin(L, a)
%
% Generate an L-point gaussian window of the given width. Use larger a
% for a narrow window.  Use larger L for a smoother curve. 
%
%     w = exp ( -(a*x)^2/2 )
%
% for x = linspace(-(L-1)/L, (L-1)/L, L)

function x = gausswin(L, w)

  if nargin < 1 || nargin > 2
    help(mfilename);
  end
  if nargin == 1, w = 2.5; end % if

  if L==1
    x = 1;
  else
    x = exp ( -0.5 * ( w/(L-1) * ( -(L-1) : 2 : L-1 )' ) .^ 2 );
  end

end
