function [nb] = find_vertex_neighbours(pnt, tri, indx)

% FIND_VERTEX_NEIGHBOURS determines the neighbours of a specified vertex
% in a mesh.
% 
% [nb] = find_vertex_neighbours(pnt, tri, indx)

% Copyright (C) 2003, Robert Oostenveld
%
% $log$

intri = any(tri==indx, 2);
nb    = tri(intri, :);
nb    = unique(nb(:));
nb    = setdiff(nb, indx);

