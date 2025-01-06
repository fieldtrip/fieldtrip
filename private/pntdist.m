function dist = pntdist(p1, p2)

% PNTDIST returns the euclidian distance between two points
%
%  [dist] = pntdist(pnt1, pnt2)
%
% where pnt1 and pnt2 must be Npnt x 3
% or either one can be Npnt x 1

% Copyright (C) 2002, Robert Oostenveld
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

if size(p1,1)==1
  p1 = repmat(p1, size(p2,1), 1);
elseif size(p2,1)==1
  p2 = repmat(p2, size(p1,1), 1);
end

dist = sqrt(sum((p1-p2).^2, 2));

