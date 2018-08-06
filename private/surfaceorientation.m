function val = surfaceorientation(pnt, tri, ori)

% SURFACEORIENTATION returns 1 if the triangulated surface is outward
% oriented, -1 if it is inward oriented and 0 if the orientation cannot be
% determined.
%
% Use as
%   surfaceorientation(pnt, tri)
% or
%   surfaceorientation(pnt, tri, ori)

% Copyright (C) 2007, Robert Oostenveld
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

if nargin<3
  ori = normals(pnt, tri, 'vertex');
end

pnt(:,1) = pnt(:,1)-mean(pnt(:,1),1);
pnt(:,2) = pnt(:,2)-mean(pnt(:,2),1);
pnt(:,3) = pnt(:,3)-mean(pnt(:,3),1);

% FIXME there is a bug in solid_angle resulting in negative values where they should be positive and vice versa 

if all(sign(sum(pnt .* ori, 2))==1)
  % the normals are outward oriented
  val = 1;
elseif all(sign(sum(pnt .* ori, 2))==-1)
  % the normals are inward oriented
  val = -1;
elseif abs(sum(solid_angle(pnt, tri))+4*pi)<1000*eps
  % the normals are outward oriented
  val = 1;
elseif abs(sum(solid_angle(pnt, tri))-4*pi)<1000*eps
  % the normals are inward oriented
  val = -1;
else
  % cannot determine
  val = 0;
end
