function [pnt, dhk] = mesh_cone(npnt)

% MESH_CONE creates a triangulated cone
%
% Use as
%   [pnt, tri] = mesh_cone(N)
%
% This creates a cone with N-2 vertices on the bottom circle and N vertices in total.
%
% See also MESH_TETRAHEDRON, MESH_OCTAHEDRON, MESH_ICOSAHEDRON, MESH_SPHERE, MESH_CUBE

% Copyright (C) 2002-2019, Robert Oostenveld
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

ncircle = npnt-2;
ndhk = 2*ncircle;

pnt = zeros(npnt, 3);
dhk = zeros(ndhk, 3);

pnt(1:ncircle, 1) = cos(linspace(0, 2*pi, ncircle)');		% bottom plane
pnt(1:ncircle, 2) = sin(linspace(0, 2*pi, ncircle)');		% bottom plane
pnt(npnt-1, :)  = [0 0 0];				% bottom point
pnt(npnt  , :)  = [0 0 1];				% top point

% connect all points on the circle to the bottom point
dhk(1:ncircle, 1) = npnt-1;
dhk(1:ncircle, 2) = (1:ncircle)';
dhk(1:ncircle, 3) = (2:(ncircle+1))';
dhk(  ncircle, 3) = 1;

% connect all points on the circle to the top point
dhk((ncircle+1):ndhk, 1) = npnt;
dhk((ncircle+1):ndhk, 2) = dhk(1:ncircle, 3);
dhk((ncircle+1):ndhk, 3) = dhk(1:ncircle, 2);

