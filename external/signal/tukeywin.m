% Copyright (C) 2007 Laurent Mazet <mazet@crm.mot.com>
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
% @deftypefn {Function File} {@var{w} =} tukeywin (@var{L}, @var{r})
% Return the filter coefficients of a Tukey window (also known as the
% cosine-tapered window) of length @var{L}. @var{r} defines the ratio
% between the constant section and and the cosine section. It has to be
% between 0 and 1. The function returns a Hanning window for @var{r}
% egals 0 and a full box for @var{r} egals 1. By default @var{r} is set
% to 1/2.
%
% For a definition of the Tukey window, see e.g. Fredric J. Harris,
% 'On the Use of Windows for Harmonic Analysis with the Discrete Fourier
% Transform, Proceedings of the IEEE', Vol. 66, No. 1, January 1978,
% Page 67, Equation 38.
% @end deftypefn

function w = tukeywin (L, r)

if nargin<2
  r = 1/2;
end

if (nargin < 1 || nargin > 2)
  help(mfilename);
elseif (nargin == 2)
  % check that 0 < r < 1
  if r > 1
    r = 1;
  elseif r < 0
    r = 0;
  end % if
end % if

% generate window
switch r
  case 0,
    % full box
    w = ones (L, 1);
  case 1,
    % Hanning window
    w = hanning (L);
  otherwise
    % cosine-tapered window
    t = linspace(0,1,L);
    t = t(1:end/2)';
    w = (1 + cos(pi*(2*t/r-1)))/2;
    w(floor(r*(L-1)/2)+2:end) = 1;
    w = [w; ones(mod(L,2)); flipud(w)];
end % switch

end

%!demo
%! L = 100;
%! r = 1/3;
%! w = tukeywin (L, r);
%! title(sprintf('%d-point Tukey window, R = %d/%d', L, [p, q] = rat(r), q));
%! plot(w);
