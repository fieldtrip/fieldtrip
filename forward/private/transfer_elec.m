function [tra] = transfer_elec(el, pos, tri)

% TRANSFER_ELEC computes the transfer matrix from vertex to electrode potentials
% using bi-linear interpolation over the triangles.
%
% Use as
%   tra = transfer_elec(el, pos, tri)
% where
%   el  = Kx3 matrix that contains the [tri, la, mu] for each electrode
%   pos = Mx3 matrix with vertex locations of the triangulated headshape
%   tri = Nx3 matrix with the vertex indices for each of the triangles
%
% See also PROJECT_ELEC, NORMALS_ELEC

% Copyright (C) 1998-2026, Robert Oostenveld
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

Npos = size(pos,1);
Ntri = size(tri,1);
Nel  = size(el,1);

tra = zeros(Nel, Npos);

for i=1:Nel
  tra(i, tri(el(i,1), 1)) = 1 - el(i,2) - el(i,3);
  tra(i, tri(el(i,1), 2)) = el(i,2);
  tra(i, tri(el(i,1), 3)) = el(i,3);
end

