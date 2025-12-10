function [side] = ptriside(v1, v2, v3, r, tolerance)

% PTRISIDE determines the side of a plane on which a set of points lie. It
% returns 0 for the points that lie exactly on the plane.
%
% [side] = ptriside(v1, v2, v3, r)
% 
% the side of points r is determined relative to the plane spanned by
% vertices v1, v2 and v3. v1,v2 and v3 should be 1x3 vectors. r should be a
% Nx3 matrix

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

if nargin<5
  tolerance = 100*eps;
end

n = size(r,1);
a = r  - ones(n,1)*v1;
b = v2 - v1;
c = v3 - v1;
d = crossproduct(b, c);
val = dotproduct(a, d);

side = zeros(n, 1);
side(val >  tolerance) =  1;
side(val < -tolerance) = -1;

%if val>tolerance
%  side=1;
%elseif val<-tolerance
%  side=-1;
%else
%  side=0;
%end

% subfunction without overhead to speed up
function c = crossproduct(a, b)

c(1) = a(2)*b(3)-a(3)*b(2);
c(2) = a(3)*b(1)-a(1)*b(3);
c(3) = a(1)*b(2)-a(2)*b(1);

% subfunction without overhead to speed up, input a can be a matrix
function d = dotproduct(a, b)

d = a(:,1)*b(1)+a(:,2)*b(2)+a(:,3)*b(3);
