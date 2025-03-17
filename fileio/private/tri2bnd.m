function [bnd, bndtissue] = tri2bnd(pos, tri, tritissue)

% TRI2BND takes a triangulated surface mesh that is represented as one
% long list of triangles with per triangle a tissue or region type, and
% converts it in a struct array with one surface mesh per tissue.
%
% Use as
%   [bnd, bndtissue] = tri2bnd(pos, tri, tritissue)
%
% The output "bnd" is a structure array with the fields bnd.pos and bnd.tri.
%
% See also MESH2EDGE, POLY2TRI, BND2TRI

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

npos = size(pos,1);
ntri = size(tri,1);

% do some sanity checks
assert(all(tri(:))>0)
assert(all(tri(:))<=npos);
assert(numel(tritissue)==ntri);

% determine the different tissue types
% note that this vector is not guaranteed to be contiguous, some tissues might be absent in the data 
bndtissue = sort(unique(tritissue));

% make a single boundary surface for each tissue type
bnd = struct();
for i=1:numel(bndtissue)
  sel = (tritissue==bndtissue(i));
  [bnd(i).pos, bnd(i).tri] = remove_unused_vertices(pos, tri(sel,:));
end
