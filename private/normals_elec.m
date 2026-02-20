function [nrm] = normals_elec(elc, pos, tri)

% NORMALS_ELEC computes the orientation of the electrodes from the given
% triangulated headshape
%
% Use as
%   nrm = normals_elec(elc, pos, tri)
% where
%   elc = Kx3 matrix with the position of all electrodes
%   pos = Mx3 matrix with vertex locations of the triangulated headshape
%   tri = Nx3 matrix with the vertex indices for each of the triangles
%
% See also PROJECT_ELEC, TRANSFER_ELEC, PCNORMALS

% Copyright (C) 2026, Robert Oostenveld
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

if strcmp(surface_orientation(pos, tri), 'inward')
  % flip the orientation of the surface
  tri = fliplr(tri);
end

nelc = size(elc,1);
npos = size(pos,1);
ntri = size(tri,1);

% to avoid confusion between electrode and headshape positions
headshape.pos = pos;
headshape.tri = tri;
clear pos tri

% compute the normals of all vertices and triangles
vertexnormals = surface_normals(headshape.pos, headshape.tri, 'vertex');
trianglenormals = surface_normals(headshape.pos, headshape.tri, 'triangle');

% determine the size of the head
headsize = prod(range(headshape.pos))^(1/3);

% take a small sphere around each electrode
radius = headsize/10;

nrm = nan(nelc,3);
for i=1:nelc
  dist = pdist2(elc(i,:), headshape.pos);
  sel = dist<radius;
  if sum(sel>10)
    nrm(i,:) = mean(vertexnormals(sel,:), 1);
    nrm(i,:) = nrm(i,:) / norm(nrm(i,:));
  else
    el = project_elec(elc(i,:), headshape.pos, headshape.tri);
    nrm(i,:) = trianglenormals(el(1),:);
  end
end

if any(isnan(nrm(:)))
  ft_warning('could not compute the orientation for each of the electrodes');
end