function [outermost, inside] = find_outermost_boundary(bnd)

% FIND_OUTERMOST_BOUNDARY locates outermost compartment of a BEM model
% by looking at the containment of the triangular meshes describing 
% the surface boundaries
%
% [outermost] = find_innermost_boundary(bnd)
%
% with the boundaries described by a struct array bnd with
%   bnd(i).pnt  vertices of boundary i (matrix of size Nx3)
%   bnd(i).tri  triangles of boundary i (matrix of size Mx3)

% Copyright (C) 2003, Robert Oostenveld

ncmp = length(bnd);

% try to locate the outermost compartment
for i=1:ncmp
for j=1:ncmp
  % determine for a single vertex on each surface if it is inside or outside the other surfaces
  curpos = bnd(i).pos(1,:); % any point on the boundary is ok
  curpnt = bnd(j).pos;
  curtri = bnd(j).tri;
  if i==j
    inside(i,j) = 0;
  else
    inside(i,j) = bounding_mesh(curpos, curpnt, curtri);
  end
end
end
% assume that the sources are in the innermost compartment
tmp = sum(inside, 2);
[i, outermost] = min(tmp);

