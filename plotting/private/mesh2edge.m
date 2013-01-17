function [newbnd] = mesh2edge(bnd)

% MESH2EDGE finds the edge lines from a triangulated mesh or the edge surfaces
% from a tetrahedral or hexahedral mesh.
%
% Use as
%   [bnd] = mesh2edge(bnd)

% Copyright (C) 2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

if isfield(bnd, 'tri')
  % make a list of all edges
  edge1 = bnd.tri(:, [1 2]);
  edge2 = bnd.tri(:, [2 3]);
  edge3 = bnd.tri(:, [3 1]);
  edge = cat(1, edge1, edge2, edge3);
  
elseif isfield(bnd, 'tet')
  % make a list of all triangles that form the tetraheder
  tri1 = bnd.tet(:, [1 2 3]);
  tri2 = bnd.tet(:, [2 3 4]);
  tri3 = bnd.tet(:, [3 4 1]);
  tri4 = bnd.tet(:, [4 1 2]);
  edge = cat(1, tri1, tri2, tri3, tri4);
  
elseif isfield(bnd, 'hex')
  % make a list of all "squares" that form the cube/hexaheder
  % FIXME should be checked, this is impossible without a drawing
  square1 = bnd.tet(:, [1 2 3 4]);
  square2 = bnd.tet(:, [5 6 7 8]);
  square3 = bnd.tet(:, [1 2 6 5]);
  square4 = bnd.tet(:, [2 3 7 6]);
  square5 = bnd.tet(:, [3 4 8 7]);
  square6 = bnd.tet(:, [4 1 5 8]);
  edge = cat(1, square1, square2, square3, square4, square5, square6);
  
end % isfield(bnd)

% make them all point in the same direction
% keep the original as "edge" and the sorted one as "sedge"
sedge = sort(edge, 2);

% find the edges that only occur once
[c, ia, ic] = unique(sedge, 'rows');
sel = false(size(ic));
for k=1:length(ic)
  sel(k) = sum(ic==k)==1;
end
% make the selection in the original, not the sorted version of the edges
% otherwise the orientation of the edges might get flipped
edge = edge(sel,:);

% the naming of the output edges depends on what they represent
newbnd.pnt  = bnd.pnt;
if isfield(bnd, 'tri')
  newbnd.edge = edge;
elseif isfield(bnd, 'tet')
  newbnd.tri = edge;
elseif isfield(bnd, 'hex')
  newbnd.poly = edge;
end
