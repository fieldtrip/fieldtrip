function [pos, tri] = remove_unused_vertices(pos, tri)

% REMOVE_UNUSED_VERTICES removes unused vertices from a triangular, tetrahedral or
% hexahedral mesh, renumbering the vertex-indices for the elements.
%
% Use as
%   [pos, tri] = remove_unused_vertices(pos, tri)
%   [pos, tet] = remove_unused_vertices(pos, tet)
%   [pos, hex] = remove_unused_vertices(pos, hex)
%
% See also REMOVE_VERTICES, REMOVE_DOUBLE_VERTICES

% Copyright (C) 2023, Robert Oostenveld
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

npos = size(pos, 1);

sel = unique(tri(:));
sel = setdiff(1:npos, sel);
[pos, tri] = remove_vertices(pos, tri, sel);
