% Copyright (C) 2008  David Bateman
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
% @deftypefn {Function File} {@var{w} =} window (@var{f}, @var{n}, @var{opts})
% Create a @var{n}-point windowing from the function @var{f}. The
% function @var{f} can be for example @code{@@blackman}. Any additional
% arguments @var{opt} are passed to the windowing function.
% @end deftypefn 

function wout = window (f, n, varargin)
  if (nargin == 0)
    error ('window: UI tool not supported');
  elseif (nargin > 1)
    w = feval (f, n, varargin{:});
    if (nargout > 0)
      wout = w;
    end % if
  else
    help(mfilename);
  end % if
end
