function val = surfaceorientation(pnt, tri, ori)

% SURFACEORIENTATION returns 1 if the triangulated surface is outward
% oriented, -1 if it is inward oriented and 0 if the orientation cannot be
% determined.
%
% Use as
%   surfaceorientation(pnt, tri)
% or
%   surfaceorientation(pnt, tri, ori)

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: surfaceorientation.m,v $
% Revision 1.1  2007/05/16 11:45:59  roboos
% new implementation
%

if nargin<3
  ori = normals(pnt, tri, 'vertex');
end

pnt(:,1) = pnt(:,1)-mean(pnt(:,1),1);
pnt(:,2) = pnt(:,2)-mean(pnt(:,2),1);
pnt(:,3) = pnt(:,3)-mean(pnt(:,3),1);

% FIXME there is a bug in solid_angle resulting in negative values where they should be positive and vice versa 

if all(sign(sum(pnt .* ori, 2))==1)
  % the normals are outward oriented
  val = 1;
elseif all(sign(sum(pnt .* ori, 2))==-1)
  % the normals are inward oriented
  val = -1;
elseif abs(sum(solid_angle(pnt, tri))+4*pi)<1000*eps
  % the normals are outward oriented
  val = 1;
elseif abs(sum(solid_angle(pnt, tri))-4*pi)<1000*eps
  % the normals are inward oriented
  val = -1;
else
  % cannot determine
  val = 0;
end
