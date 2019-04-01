function [pos, tri] = icosahedron(varargin)

% ICOSAHEDRON creates an icosahedron with 12 vertices and 20 triangles
%
% Use as
%   [pos, tri] = icosahedron
%
% See also TETRAHEDRON, OCTAHEDRON, ICOSAHEDRON42, ICOSAHEDRON162, ICOSAHEDRON642, ICOSAHEDRON2562

% Copyright (C) 2002, Robert Oostenveld
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

warning('icosahedron is only a compatibility wrapper, which will soon be removed. Please instead call sphere_mesh.');
[pos, tri] = sphere_mesh(varargin{:});
