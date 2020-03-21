function [pos, tri] = mesh_octahedron

% MESH_OCTAHEDRON returns the vertices and triangles of an octahedron
%
% Use as
%   [pos tri] = mesh_octahedron;
%
% See also MESH_TETRAHEDRON, MESH_OCTAHEDRON, MESH_SPHERE

% Copyright (C) 2019, Robert Oostenveld and Jan-Mathijs Schoffelen
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

pos = [
  0  0 +1
  +1  0  0
  0 +1  0
  -1  0  0
  0 -1  0
  0  0 -1
  ];

tri = [
  1 2 3
  1 3 4
  1 4 5
  1 5 2
  6 3 2
  6 4 3
  6 5 4
  6 2 5
  ];

  % scale all vertices to the unit sphere
  pos = pos ./ repmat(sqrt(sum(pos.^2,2)), 1,3);
  
