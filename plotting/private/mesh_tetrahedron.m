function [pos, tri] = mesh_tetrahedron(n)

% MESH_TETRAHEDRON returns the vertices and triangles of a refined tetrahedron,
% with n refinement steps.
%
% Use as
%   [pos, tri] = mesh_tetrahedron(N);
%
% See also MESH_ICOSAHEDRON, MESH_OCTAHEDRON, MESH_SPHERE

% Copyright (C) 2018-2019, Robert Oostenveld and Jan-Mathijs Schoffelen
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

if nargin==0 || isempty(n)
  n = 0;
end

v1 = [ 0, 0, 1 ];
v2 = [  sqrt(8/9),          0, -1/3 ];
v3 = [ -sqrt(2/9),  sqrt(2/3), -1/3 ];
v4 = [ -sqrt(2/9), -sqrt(2/3), -1/3 ];

pos = [
  v1
  v2
  v3
  v4
  ];

tri = [
  1 2 3
  1 3 4
  1 4 2
  2 4 3
  ];

if n>0
  % perform an n-fold refinement
  for i=1:n
    [pos, tri] = refine(pos, tri, 'banks');
  end
  % scale all vertices to the unit sphere
  pos = pos ./ repmat(sqrt(sum(pos.^2,2)), 1,3);
end
