function [area] = surface_area(pos, tri)

% SURFACE_AREA computes the surface area of each of the triangles in a mesh
%
% Use as
%   area = surface_area(pos, tri)
%
% See also SURFACE_ORIENTATION, SURFACE_INSIDE, SURFACE_NESTING, PROJECTTRI, PCNORMALS

% Copyright (C) 2024, Robert Oostenveld
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

% see https://www.maths.usyd.edu.au/u/MOW/vectors/vectors-11/v-11-7.html
vec12 = pos(tri(:,2),:)-pos(tri(:,1),:);
vec13 = pos(tri(:,3),:)-pos(tri(:,1),:);
area = sqrt(sum(cross(vec12, vec13).^2, 2))/2;
