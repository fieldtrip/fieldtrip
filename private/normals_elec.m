function [nrm] = normals_elec(elc, pos, tri)

% NORMALS_ELEC computes the surface normals for the electrodes from the given
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

headshape.pos = pos;
headshape.tri = tri;
headshape.nrm = surface_normals(headshape.pos, headshape.tri);

% compute the orientation for each electrode
dist = pdist2(elc, elc);
mindist = zeros(1, size(elc,1));
for i=1:size(elc,1)
  d = dist(i,:);
  d(i) = [];
  mindist(i) = min(d);
end

% determine the size of the head
headsize = prod(range(headshape.pos))^(1/3);

% take a small sphere around each electrode
radius = headsize/10;

nrm = zeros(size(elc,1),3);
for i=1:size(elc,1)
  dist = pdist2(elc(i,:), headshape.pos);
  sel = dist<radius;
  nrm(i,:) = mean(headshape.nrm(sel,:), 1);
end
