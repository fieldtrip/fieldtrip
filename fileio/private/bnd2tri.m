function [pos, tri, tissue] = bnd2tri(bnd)

% BND2TRI takes a struct array with one triangulated surface mesh per
% tissue and converts it into a single triangulated surface mesh
% represented as one long list of triangles with per triangle a tissue or
% region type.
%
% Use as
%   [pos, tri, tissue] = bnd2tri(bnd)
%
% See also MESH2EDGE, POLY2TRI, TRI2BND

% Copyright (C) 2025, Robert Oostenveld
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


pos = zeros(0,3);
tri = zeros(0,3);
tissue = zeros(0,1);

for i=1:numel(bnd)
  ntri = size(bnd(i).tri);  % of this boundary
  offset = size(pos,1);     % of all vertices together
  pos = cat(1, pos, bnd(i).pos);
  tri = cat(1, pos, bnd(i).tri + offset);
  tissue = cat(1, tissue, i*ones(ntri,1));
end

[newpos, newtri] = remove_double_vertices(pos, tri);

if ~isequal(tri, newtri)
  % this happens when the vertex sorting is not stable
  clear tissue
end

pos = newpos;
tri = newtri;
