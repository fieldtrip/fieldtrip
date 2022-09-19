function [pnt, tri] = mesh_cube()

% MESH_CUBE creates a triangulated cube
%
% Use as
%   [pos, tri] = mesh_cube()
%
% See also MESH_TETRAHEDRON, MESH_OCTAHEDRON, MESH_ICOSAHEDRON, MESH_SPHERE, MESH_CONE

% Copyright (C) 2019, Robert Oostenveld
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

pnt = [
  -1 -1 -1
  -1  1 -1
   1  1 -1
   1 -1 -1
  -1 -1  1
  -1  1  1
   1  1  1
   1 -1  1
  ];

tri = [
  1 2 4
  2 3 4
  1 5 6
  1 6 2
  2 6 7
  2 7 3
  3 7 8
  3 8 4
  4 8 5
  4 5 1
  8 6 5
  8 7 6
  ];
