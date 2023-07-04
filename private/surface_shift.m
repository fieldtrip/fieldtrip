function [pos] = surface_shift(pos, tri, amount)

% SURFACE_SHIFT inflates or deflates a triangulated surface by moving the
% vertices outward or inward along their normals.
%
% Use as
%   pos = surface_inflate(pos, tri, amount)
% where pos and tri describe the surface.
%
% See also SURFACE_NORMALS, SURFACE_ORIENTATION, SURFACE_INSIDE,
% SURFACE_NESTING

if isempty(tri)
  propos = elproj(pos);   % projection to 2D
  tri = delaunay(propos); % creates delaunay triangulation of 2D plane, which will be used for the the 3D case
end

nrm = surface_normals(pos, tri); % compute normals of surface
switch surface_orientation(pos, tri, nrm)
  case 'outward'
    % move along the normals
    nrm = +nrm;
  case 'inward'
    % move opposite to the normals
    nrm = -nrm;
  otherwise
    ft_warning('cannot determine the orientation of the surface');
end

% move the vertices with the specified amount along the normals
pos = pos + amount*nrm;
