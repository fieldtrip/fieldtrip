function [inside] = bounding_mesh(pos, pnt, tri);

% BOUNDING_MESH determines if a point is inside/outside a triangle mesh 
% whereby the bounding triangle mesh should be closed.
%
% [inside] = bounding_mesh(pos, pnt, tri)
%
% where
%   pos     position of point of interest (can be 1x3 or Nx3)
%   pnt     bounding mesh vertices
%   tri     bounding mesh triangles
%
% See also SOLID_ANGLE

% Copyright (C) 2003, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
% $Id: bounding_mesh.m 2885 2011-02-16 09:41:58Z roboos $

global fb;
if isempty(fb)
  fb = 0;
end

npos = size(pos, 1);
npnt = size(pnt, 1);
ntri = size(tri, 1);

% determine a cube that encompases the boundary triangulation
bound_min = min(pnt);
bound_max = max(pnt);

% determine a sphere that is completely inside the boundary triangulation
bound_org = mean(pnt);
bound_rad = sqrt(min(sum((pnt - repmat(bound_org, size(pnt,1), 1)).^2, 2)));

inside = zeros(npos, 1);
for i=1:npos
  if fb
    fprintf('%6.2f%%', 100*i/npos);
  end
  if any(pos(i,:)<bound_min) || any(pos(i,:)>bound_max)
    % the point is outside the bounding cube
    inside(i) = 0;
    if fb, fprintf(' outside the bounding cube\n'); end
  elseif sqrt(sum((pos(i,:)-bound_org).^2, 2))<bound_rad 
    % the point is inside the interior sphere
    inside(i) = 1;
    if fb, fprintf(' inside the interior sphere\n'); end
  else
    % the point is inside the bounding cube but outside the interior sphere
    % compute the total solid angle of the surface, which is zero for a point outside
    % the triangulation and 4*pi or -4*pi for a point inside (depending on the triangle
    % orientation)
    tmp = pnt - repmat(pos(i,:), npnt, 1);
    solang = solid_angle(tmp, tri);
    if any(isnan(solang))
      inside(i) = nan;
    elseif (abs(sum(solang))-2*pi)<0
      % total solid angle is (approximately) zero
      inside(i) = 0;
    elseif (abs(sum(solang))-2*pi)>0
      % total solid angle is (approximately) plus or minus 4*pi
      inside(i) = 1;
    end
    if fb, fprintf(' solid angle\n'); end
  end
end

