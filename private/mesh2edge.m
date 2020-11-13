function [output] = mesh2edge(mesh)

% MESH2EDGE finds the edge lines from a triangulated mesh or the edge
% surfaces from a tetrahedral or hexahedral mesh. An edge is defined as an
% element that does not border any other element. This also implies that a
% closed triangulated surface has no edges.
%
% Use as
%   [edge] = mesh2edge(mesh)
%
% See also POLY2TRI

% Copyright (C) 2013-2020, Robert Oostenveld
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

if isfield(mesh, 'tri')
  % make a list of all edges
  edge1 = mesh.tri(:, [1 2]);
  edge2 = mesh.tri(:, [2 3]);
  edge3 = mesh.tri(:, [3 1]);
  edge = cat(1, edge1, edge2, edge3);
  
elseif isfield(mesh, 'tet')
  % make a list of all triangles that form the tetraheder
  tri1 = mesh.tet(:, [1 2 3]);
  tri2 = mesh.tet(:, [2 3 4]);
  tri3 = mesh.tet(:, [3 4 1]);
  tri4 = mesh.tet(:, [4 1 2]);
  edge = cat(1, tri1, tri2, tri3, tri4);
  
elseif isfield(mesh, 'hex')
  % make a list of all "squares" that form the cube/hexaheder
  % FIXME should be checked, this is impossible without a drawing
  square1 = mesh.hex(:, [1 2 3 4]);
  square2 = mesh.hex(:, [5 6 7 8]);
  square3 = mesh.hex(:, [1 2 6 5]);
  square4 = mesh.hex(:, [2 3 7 6]);
  square5 = mesh.hex(:, [3 4 8 7]);
  square6 = mesh.hex(:, [4 1 5 8]);
  edge = cat(1, square1, square2, square3, square4, square5, square6);
  
end % isfield(mesh)

% soort all polygons in the same direction
% keep the original as "edge" and the sorted one as "sedge"
sedge = sort(edge, 2);

% % find the edges that are not shared -> count the number of occurences
% n = size(sedge,1);
% occurences = ones(n,1);
% for i=1:n
%   for j=(i+1):n
%     if all(sedge(i,:)==sedge(j,:))
%       occurences(i) = occurences(i)+1;
%       occurences(j) = occurences(j)+1;
%     end
%   end
% end
%
% % make the selection in the original, not the sorted version of the edges
% % otherwise the orientation of the edges might get flipped
% edge = edge(occurences==1,:);

% find the edges that are not shared
indx = findsingleoccurringrows(sedge);
edge = edge(indx, :);

% replace pnt by pos
mesh = fixpos(mesh);

% the naming of the edges in the output depends on what they represent
output.pos = mesh.pos;
if isfield(mesh, 'tri')
  % these have two vertices in each edge element
  output.line = edge;
elseif isfield(mesh, 'tet')
  % these have three vertices in each edge element
  output.tri = edge;
elseif isfield(mesh, 'hex')
  % these have four vertices in each edge element
  output.poly = edge;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION, see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1833#c12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function indx = findsingleoccurringrows(X)
[X, indx] = sortrows(X);
sel  = any(diff([X(1,:)-1; X],1),2) & any(diff([X; X(end,:)+1],1),2);
indx = indx(sel);

function indx = finduniquerows(X)
[X, indx] = sortrows(X);
sel  = any(diff([X(1,:)-1; X],1),2);
indx = indx(sel);
