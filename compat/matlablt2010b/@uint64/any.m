function y = any(x)

% ANY True if any element of a vector is a nonzero number or is
% logical 1 (TRUE).  ANY ignores entries that are NaN (Not a Number).
%
% For vectors, ANY(V) returns logical 1 (TRUE) if any of the
% elements of the vector is a nonzero number or is logical 1 (TRUE).
% Otherwise it returns logical 0 (FALSE).  For matrices, ANY(X)
% operates on the columns of X, returning a row vector of logical 1's
% and 0's.

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin>1
  error('this implementation is only supported with one input argument');
end

siz = size(x);
if numel(siz)>2
  error('this implementation is only supported with vector or matrix input');
end

if siz(1)==1 || siz(2)==1
  y = false;
  for i=1:prod(siz)
    if x(i)
      y = true;
      break
    end
  end

else
  y = false(1,siz(2));
  for j=1:siz(2)
    for i=1:siz(1)
      if x(i,j)
        y(1,j) = true;
        break
      end
    end
  end

end

