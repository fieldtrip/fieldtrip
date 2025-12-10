function s = surface_orientation(pos, tri, ori)

% SURFACE_ORIENTATION returns the string 'inward' or 'outward' or 'unknown',
% depending on the surface orientation.
%
% Use as
%   str = surface_orientation(pos, tri)
% or
%   str = surface_orientation(pos, tri, ori)
%
% See also SURFACE_AREA, SURFACE_NESTING, SURFACE_NORMALS, SURFACE_NESTING

% Copyright (C) 2007-2024, Robert Oostenveld
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

assert(isnumeric(pos));
assert(isnumeric(tri));

% translate the center to the origin
center = mean(pos,1);
pos(:,1) = pos(:,1) - center(1);
pos(:,2) = pos(:,2) - center(2);
pos(:,3) = pos(:,3) - center(3);

if nargin==3
  % look at the orientation of the normals seen from the center
  % this method is rigorous only for star-shaped surfaces, but it does also work for open surfaces
  n = sign(sum(pos .* ori, 2));
  
  if all(n==1)
    s = 'outward';
    return
  elseif all(n==-1)
    s = 'inward';
    return
  end
end % if normals provided

% look at the sum of the solid angles as seen from the center
% it assumes that the center of the surface is inside
% this is computationally more expensive and requires the surface to be closed
solang = sum(solid_angle(pos, tri));

tolerance = 1000*eps;

if solang<0 && (abs(solang)-4*pi)<tolerance
  s = 'outward';
elseif solang>0 && (abs(solang)-4*pi)<tolerance
  s = 'inward';
else
  s = 'unknown';
end
