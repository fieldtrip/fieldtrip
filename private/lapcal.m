function [lap, edge] = lapcal(pnt, tri)

% LAPCAL computes the finite difference approximation to the surface laplacian
% matrix using a triangulation of the surface
%
% lap = lapcal(pnt, tri)
%
% where
%   pnt	  contains the positions of the vertices
%   tri	  contains the triangle definition
%   lap	  is the surface laplacian matrix
%
% See also LAPINT, LAPINTMAT, READ_TRI, SAVE_TRI

% For details see
%   T.F. Oostendorp, A. van Oosterom, and G.J.M. Huiskamp.
%   Interpolation on a triangulated 3D surface.
%   Journal of Computational Physics, 80:331-343, 1989. 

% Copyright (C) 2001, Robert Oostenveld
%
% $Log: lapcal.m,v $
% Revision 1.1  2009/01/21 10:46:10  roboos
% moved from forwinv/* and forwinv/mex/* directory to forwinv/private/* to make the CVS layout consistent with the release version
%
% Revision 1.2  2003/03/11 14:45:36  roberto
% updated help and copyrights
%

npnt = size(pnt,1);
ntri = size(tri,1);

% the matrix edge is the connectivity matric for all vertices
edge = spalloc(npnt, npnt, 3*ntri);
for i=1:ntri
  % compute the lenght of all triangle edges
  edge(tri(i,1), tri(i,2)) = norm(pnt(tri(i,1),:) - pnt(tri(i,2),:));
  edge(tri(i,2), tri(i,3)) = norm(pnt(tri(i,2),:) - pnt(tri(i,3),:));
  edge(tri(i,3), tri(i,1)) = norm(pnt(tri(i,3),:) - pnt(tri(i,1),:));
  % make sure that all edges are symmetric
  edge(tri(i,2), tri(i,1)) = edge(tri(i,1), tri(i,2));
  edge(tri(i,3), tri(i,2)) = edge(tri(i,2), tri(i,3));
  edge(tri(i,1), tri(i,3)) = edge(tri(i,3), tri(i,1));
end

lap = zeros(npnt);
for i=1:npnt
  k = find(edge(i,:));		% the indices of the neighbours
  n = length(k);		% the number of neighbours
  h  = mean(edge(i,k));		% the average distance to the neighbours
  hi = mean(1./edge(i,k));	% the average inverse distance to the neighbours

  lap(i,i) = -(4/h) * hi;
  lap(i,k) =  (4/(h*n)) * 1./edge(i,k);
end

