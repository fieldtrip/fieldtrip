function [pnt, line] = find_mesh_edge(pnt, tri)

% FIND_MESH_EDGE returns the edge of a triangulated mesh
%
% [pnt, line] = find_mesh_edge(pnt, tri), where
%
% pnt   contains the vertex locations and 
% line  contains the indices of the linepieces connecting the vertices

% Copyright (C) 2003, Robert Oostenveld
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

npnt = size(pnt,1);
ntri = size(tri,1);

line = zeros(0,2);
nb   = find_triangle_neighbours(pnt, tri);
edge = isnan(nb);
indx = any(edge, 2);
for i=find(indx)'
  if isnan(nb(i,1))
    line(end+1, :) = tri(i,[1 2]);
  end
  if isnan(nb(i,2))
    line(end+1, :) = tri(i,[2 3]);
  end
  if isnan(nb(i,3))
    line(end+1, :) = tri(i,[3 1]);
  end
end
