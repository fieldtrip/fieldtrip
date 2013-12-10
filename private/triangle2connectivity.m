function [connmat] = triangle2connectivity(tri)

% TRIANGLE2CONNECTIVITY computes a connectivity-matrix from a triangulation.
%
% Use as
%  [connmat] = triangle2connectivity(tri)
%
% The input tri is an Nx3 matrix describing a triangulated surface,
% containing indices to connecting vertices
% The output connmat is a sparse logical NxN matrix, with ones, where vertices
% are connected, and zeros otherwise.

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
connmat  = sparse([edges(:,1);edges(:,2)],[edges(:,2);edges(:,1)],true(2*n,1));

function [newtri] = tri_reindex(tri)

% this subfunction reindexes tri such that they run from 1:number of unique vertices
newtri       = tri;
[srt, indx]  = sort(tri(:));
tmp          = cumsum(double(diff([0;srt])>0));
newtri(indx) = tmp;
