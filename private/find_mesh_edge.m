function [pnt, line] = find_mesh_edge(pnt, dhk);

% FIND_MESH_EDGE returns the edge of a triangulated mesh
%
% [pnt, line] = find_mesh_edge(pnt, dhk), where
%
% pnt	contains the vertex locations and 
% line	contains the indices of the linepieces connecting the vertices

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: find_mesh_edge.m,v $
% Revision 1.2  2003/03/04 21:46:18  roberto
% added CVS log entry and synchronized all copyright labels
%

npnt = size(pnt,1);
ndhk = size(dhk,1);

line = zeros(0,2);
nb   = find_triangle_neighbours(pnt, dhk);
edge = isnan(nb);
indx = any(edge, 2);
for i=find(indx)'
  if isnan(nb(i,1))
    line(end+1, :) = dhk(i,[1 2]);
  end
  if isnan(nb(i,2))
    line(end+1, :) = dhk(i,[2 3]);
  end
  if isnan(nb(i,3))
    line(end+1, :) = dhk(i,[3 1]);
  end
end
