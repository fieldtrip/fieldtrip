function [posR, triR, removetri] = remove_vertices(pos, tri, removepos)

% REMOVE_VERTICES removes specified indexed vertices from a triangular mesh
% renumbering the vertex-indices for the triangles and removing all
% resulting 'open' triangles.
%
% Use as
%   [pos, tri] = remove_vertices(pos, tri, sel)

% Copyright (C) 2004-2022, Robert Oostenveld
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

if numel(removepos)==size(pos,1) && all(removepos==0 | removepos==1)
  removepos = find(removepos);
end

% remove the vertices and determine the new numbering (indices) in numb
keeppos = setdiff(1:npos, removepos);
numb    = zeros(1,npos);
numb(keeppos) = 1:length(keeppos);

% look for triangles referring to removed vertices
removetri = false(ntri,1);
for i=1:size(tri,2)
  % loop over all columns of tri, this also works for tet and hex
  removetri(ismember(tri(:,i), removepos)) = true;
end

% remove the vertices and triangles
posR = pos(keeppos, :);
triR = tri(~removetri,:);

% renumber the vertex indices for the triangles
triR = numb(triR);
