% Copyright (C) 2000-2002 Paul Kienzle <pkienzle@users.sf.net>
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
% usage:  w = triang (L)
%
% Returns the filter coefficients of a triangular window of length L.
% Unlike the bartlett window, triang does not go to zero at the edges
% of the window.  For odd L, triang(L) is equal to bartlett(L+2) except
% for the zeros at the edges of the window.

function w = triang(L)
  if (nargin ~= 1)
    help(mfilename);
  elseif (~isscalar(L) || L ~= fix (L) || L < 1)
    error('triang: L has to be an integer > 0');
  end % if
  w = 1 - abs ([-(L-1):2:(L-1)]' / (L+rem(L,2)));
end

%!error triang
%!error triang(1,2)
%!error triang([1,2]);
%!assert (triang(1), 1)
%!assert (triang(2), [1; 1]/2)
%!assert (triang(3), [1; 2; 1]/2);
%!assert (triang(4), [1; 3; 3; 1]/4);
%!test
%! x = bartlett(5);
%! assert (triang(3), x(2:4));

%!demo
%! subplot(221); axis([-1, 1, 0, 1.3]); grid('on');
%! title('comparison with continuous for odd n');
%! n=7; k=(n-1)/2; t=[-k:0.1:k]/(k+1);
%! plot(t,1-abs(t),';continuous;',[-k:k]/(k+1),triang(n),'g*;discrete;');
%!
%! subplot(222); axis([-1, 1, 0, 1.3]); grid('on');
%! n=8; k=(n-1)/2; t=[-k:0.1:k]/(k+1/2);
%! title('note the higher peak for even n');
%! plot(t,1+1/n-abs(t),';continuous;',[-k:k]/(k+1/2),triang(n),'g*;discrete;');
%!
%! subplot(223); axis; grid('off');
%! title('n odd, triang(n)==bartlett(n+2)');
%! n=7;
%! plot(0:n+1,bartlett(n+2),'g-*;bartlett;',triang(n),'r-+;triang;');
%!
%! subplot(224); axis; grid('off');
%! title('n even, triang(n)!=bartlett(n+2)');
%! n=8;
%! plot(0:n+1,bartlett(n+2),'g-*;bartlett;',triang(n),'r-+;triang;');
%!
%! subplot(111); title('');
