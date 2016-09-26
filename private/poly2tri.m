function mesh = poly2tri(mesh)

% POLY2TRI converts the polygons in a mesh to triangles by splitting
% them in half. The input polygons should consist of 4 vertices.
% Curvature is not considered and the resulting split will only be
% optimal for flat polygons.
%
% Use as
%  mesh = poly2tri(mesh)
%
% See also MESH2EDGE

if ~isfield(mesh, 'poly')
  error('the mesh does not contain polygons')
elseif ~size(mesh.poly,2)==4
  error('this is only implemented for polygons with 4 corner points')
end

% there are 4 vertices per polygon, which do not have to ly in a single plane
% hence there are two possilbe triangulations to split the polygon [1 2 3 4]
% option 1 is [1 2 3] and [1 3 4] -> this is used in the following code
% option 2 is [4 1 2] and [4 2 3]

npoly = size(mesh.poly,1);

% make a vectorized representation of the polygons (i.e. vertex indices)
poly = reshape(mesh.poly', npoly*4, 1);

% determine the vertex indices for the triangles
tri = [
  (1:4:npoly*4)' (2:4:npoly*4)' (3:4:npoly*4)'
  (1:4:npoly*4)' (3:4:npoly*4)' (4:4:npoly*4)'
  ];

% look up the triangle indices from the polygon indices
mesh.tri = poly(tri);

% remove the polygons
mesh = rmfield(mesh, 'poly');
