function [pos, tri, tissue] = remove_double_vertices(pos, tri, tissue)

% REMOVE_DOUBLE_VERTICES removes double vertices from a triangular, tetrahedral or
% hexahedral mesh, renumbering the vertex-indices for the elements.
%
% Use as
%   [pos, tri] = remove_double_vertices(pos, tri)
%   [pos, tet] = remove_double_vertices(pos, tet)
%   [pos, hex] = remove_double_vertices(pos, hex)
%
% See also REMOVE_VERTICES, REMOVE_UNUSED_VERTICES

% Copyright (C) 2004-2025, Robert Oostenveld and Jan-Mathijs Schoffelen
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

npos = size(pos,1);
ntri = size(tri,1);

% do some sanity checks
assert(all(tri(:))>0)
assert(all(tri(:))<=npos);
if nargin>2
  assert(numel(tissue)==ntri);
end

% find unique vertices, keep them in the same order
[newpos, ia, ic] = unique(pos, 'rows', 'stable');

% determine the mapping from the old to new indices
mapping = zeros(1,npos);
mapping(ia) = 1:length(ia);

% re-index the vertex indices that are represented in tri
tri = ia(ic(tri));
tri = mapping(tri);

% return only the unique vertices
pos = newpos;
