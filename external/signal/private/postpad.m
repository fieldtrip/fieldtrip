% Copyright (C) 1994, 1995, 1996, 1997, 1998, 2000, 2002, 2004, 2005,
%               2006, 2007, 2008, 2009 John W. Eaton
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
% @deftypefn {Function File} {} postpad (@var{x}, @var{l}, @var{c})
% @deftypefnx {Function File} {} postpad (@var{x}, @var{l}, @var{c}, @var{dim})
% @seealso{prepad, resize}
% @end deftypefn

% Author: Tony Richardson <arichard@stark.cc.oh.us>
% Created: June 1994

function y = postpad (x, l, c, dim)

if nargin < 2 || nargin > 4
  %print_usage ();
  error('wrong number of input arguments, should be between 2 and 4');
end

if nargin < 3 || isempty(c)
  c = 0;
else
  if ~isscalar(c)
    error ('postpad: third argument must be empty or a scalar');
  end
end

nd = ndims(x);
sz = size(x);
if nargin < 4
  % Find the first non-singleton dimension
  dim  = 1;
  while dim < nd+1 && sz(dim)==1
    dim = dim + 1;
  end
  if dim > nd
    dim = 1;
  elseif ~(isscalar(dim) && dim == round(dim)) && dim > 0 && dim< nd+1
    error('postpad: dim must be an integer and valid dimension');
  end
end

if ~isscalar(l) || l<0
  error ('second argument must be a positive scalar');
end

if dim > nd
  sz(nd+1:dim) = 1;
end

d = sz(dim);

if d >= l
  idx = cell(1,nd);
  for i = 1:nd
    idx{i} = 1:sz(i);
  end
  idx{dim} = 1:l;
  y = x(idx{:});
else
  sz(dim) = l-d;
  y = cat(dim, x, c * ones(sz));
end

