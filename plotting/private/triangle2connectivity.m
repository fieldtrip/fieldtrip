function [connmat] = triangle2connectivity(tri, pos)

% TRIANGLE2CONNECTIVITY computes a connectivity-matrix from a triangulation.
%
% Use as
%  [connmat] = triangle2connectivity(tri)
% or
%  [connmat] = triangle2connectivity(tri, pos)
%
% The input tri is an Mx3 matrix describing a triangulated surface,
% containing indices to connecting vertices. The output connmat is a sparse
% logical NxN matrix, with ones, where vertices are connected, and zeros
% otherwise.
%
% If you specify the vertex positions in the second input argument as Nx3
% matrix, the output will be a sparse matrix with the lengths of the
% edges between the connected vertices.

% Copyright (C) 2015, Jan-Mathijs Schoffelen
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

% ensure that the vertices are indexed starting from 1
if min(tri(:))==0,
  tri = tri + 1;
end

% ensure that the vertices are indexed according to 1:number of unique vertices
tri = tri_reindex(tri);

% create the unique edges from the triangulation
edges  = [tri(:,[1 2]); tri(:,[1 3]); tri(:,[2 3])];
edges  = double(unique(sort([edges; edges(:,[2 1])],2), 'rows'));

% fill the connectivity matrix
n        = size(edges,1);
if nargin<2
  % create sparse binary matrix
  connmat = sparse([edges(:,1);edges(:,2)],[edges(:,2);edges(:,1)],true(2*n,1));
else
  % create sparse matrix with edge lengths
  dpos    = sqrt(sum( (pos(edges(:,1),:) - pos(edges(:,2),:)).^2, 2));
  connmat = sparse([edges(:,1);edges(:,2)],[edges(:,2);edges(:,1)],[dpos(:);dpos(:)]);
end

function [newtri] = tri_reindex(tri)

% this subfunction reindexes tri such that they run from 1:number of unique vertices
newtri       = tri;
[srt, indx]  = sort(tri(:));
tmp          = cumsum(double(diff([0;srt])>0));
newtri(indx) = tmp;
