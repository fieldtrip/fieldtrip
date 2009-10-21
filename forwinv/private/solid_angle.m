function [w] = solid_angle(r1, r2, r3);

% SOLID_ANGLE of a planar triangle as seen from the origin
%
% The solid angle W subtended by a surface S is defined as the surface
% area W of a unit sphere covered by the surface's projection onto the
% sphere. Solid angle is measured in steradians, and the solid angle
% corresponding to all of space being subtended is 4*pi sterradians.
% 
% Use: 
%   [w] = solid_angle(v1, v2, v3) or
%   [w] = solid_angle(pnt, tri)
% where v1, v2 and v3 are the vertices of a single triangle in 3D or
% pnt and tri contain a description of a triangular mesh (this will
% compute the solid angle for each triangle)
% 
% also implemented as MEX file

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: solid_angle.m,v $
% Revision 1.3  2003/03/21 13:32:53  roberto
% created mex implementation
%
% Revision 1.2  2003/03/04 21:46:19  roberto
% added CVS log entry and synchronized all copyright labels
%

if nargin==2
  % reassign the input arguments 
  pnt = r1;
  tri = r2;
  npnt = size(pnt,1);
  ntri = size(tri,1);
  w    = zeros(ntri,1);
  % compute solid angle for each triangle
  for i=1:ntri
    r1 = pnt(tri(i,1),:);
    r2 = pnt(tri(i,2),:);
    r3 = pnt(tri(i,3),:);
    w(i) = solid_angle(r1, r2, r3);
  end 
  return
elseif nargin==3
  % compute the solid angle for this triangle  
  cp23_x = r2(2) * r3(3) - r2(3) * r3(2);
  cp23_y = r2(3) * r3(1) - r2(1) * r3(3);
  cp23_z = r2(1) * r3(2) - r2(2) * r3(1);
  nom = cp23_x * r1(1) + cp23_y * r1(2) + cp23_z * r1(3);
  n1 = sqrt (r1(1) * r1(1) + r1(2) * r1(2) + r1(3) * r1(3));
  n2 = sqrt (r2(1) * r2(1) + r2(2) * r2(2) + r2(3) * r2(3));
  n3 = sqrt (r3(1) * r3(1) + r3(2) * r3(2) + r3(3) * r3(3));
  ip12 = r1(1) * r2(1) + r1(2) * r2(2) + r1(3) * r2(3);
  ip23 = r2(1) * r3(1) + r2(2) * r3(2) + r2(3) * r3(3);
  ip13 = r1(1) * r3(1) + r1(2) * r3(2) + r1(3) * r3(3);
  den = n1 * n2 * n3 + ip12 * n3 + ip23 * n1 + ip13 * n2; 
  if (nom == 0)
    if (den <= 0)
      w = nan;
      return
    end
  end
  w = 2 * atan2 (nom, den);
  return
else
  error('invalid input');
end

