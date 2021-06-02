 function [X, Y, Z] = extract_contour(pos, tri, indx, mask)

% EXTRACT_BOUNDARY generates contour around specified vertices of a
% triangulated surface.
% It returns the coordinates of the vertices which form the contour

% % Use as
%   [X, Y, Z] = extract_contour(pos, tri, indx, mask)
%
% where mask specifies which vertices are inside or outside the contour
% and indx contains indices of which vertices should be group within a
% single contour. The return values are each Nx2 for the N line segments making up the contour.

% Copyright (C) Sophie Arana
%
% $Id$

outbnd = [];
for i = 1:length(indx)
    % for each vertex in mask
    [row,dum] = find(tri == (indx(i))); % find neighbours
    neigh = tri(row,:);
    neigh = unique(neigh(:));
    outbnd = [outbnd neigh(mask(neigh)==0)'];
    outbnd(outbnd==indx(i)) = [];
end % find all "outer" neighbours to this cluster
outbnd = unique(outbnd);

indx1 = any(ismember(tri,outbnd)'); % outer vertex connections
indx2 = any(ismember(tri,indx)'); % inner vertex connections
indxnew = find(indx1 & indx2); % triangles spanning at least one outer & one inner vertex

newpos = [];
for  i = 1:length(indxnew)
    idx = ismember(tri(indxnew(i),:),indx);
    newpos = [newpos ; pos(tri(indxnew(i),idx),:) - ((pos(tri(indxnew(i),idx),:) - pos(tri(indxnew(i),~idx),:))*0.5)];
end
n = size(newpos,1);
X = reshape(newpos(:,1),[2 n/2])';
Y = reshape(newpos(:,2),[2 n/2])';
Z = reshape(newpos(:,3),[2 n/2])';