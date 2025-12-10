function [inside] = surface_inside(dippos, pos, tri)

% SURFACE_INSIDE determines if a point is inside/outside a triangle mesh
% whereby the bounding triangle mesh should be closed.
%
% Use as
%   inside = surface_inside(dippos, pos, tri)
% where
%   dippos  position of point of interest (can be 1x3 or Nx3)
%   pos     bounding mesh vertices
%   tri     bounding mesh triangles
%
% See also SURFACE_AREA, SURFACE_ORIENTATION, SURFACE_NORMALS, SURFACE_NESTING, SOLID_ANGLE

% Copyright (C) 2003, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistrinside = surface_inside(dippos, pos, tri)ibute it and/or modify
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

% this is a work-around for http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2369
dippos = double(dippos);
pos = double(pos);
tri = double(tri);

ndip = size(dippos, 1);
npos = size(pos, 1);
ntri = size(tri, 1);

% determine a cube that encompases the boundary triangulation
cube_min = min(pos);
cube_max = max(pos);

% determine a sphere that is completely inside the boundary triangulation
sphere_center = mean(pos);
sphere_radius = sqrt(min(sum((pos - repmat(sphere_center, size(pos,1), 1)).^2, 2)));

tolerance = 1000*eps;
inside = zeros(ndip, 1);

for i=1:ndip
  if any(dippos(i,:)<cube_min) || any(dippos(i,:)>cube_max)
    % the point is outside the bounding cube
    inside(i) = 0;
  elseif sqrt(sum((dippos(i,:)-sphere_center).^2, 2))<sphere_radius
    % the point is inside the interior sphere
    inside(i) = 1;
  else
    % the point is inside the bounding cube but outside the interior sphere
    % compute the total solid angle of the surface, which is zero for a point outside
    % the triangulation and 4*pi or -4*pi for a point inside (depending on the triangle
    % orientation)
    tmp = pos - repmat(dippos(i,:), npos, 1);
    solang = solid_angle(tmp, tri);
    if any(isnan(solang))
      inside(i) = nan;
    elseif abs(sum(solang)) < tolerance
      % total solid angle is (approximately) zero
      inside(i) = 0;
    elseif (abs(sum(solang))-4*pi) < tolerance
      % total solid angle is (approximately) plus or minus 4*pi
      inside(i) = 1;
    end
  end
end
