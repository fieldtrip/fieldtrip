%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = surface_orientation(bnd)
points = bnd.pos;
faces = bnd.tri;

% translate to the center
org = mean(points,1);
points(:,1) = points(:,1) - org(1);
points(:,2) = points(:,2) - org(2);
points(:,3) = points(:,3) - org(3);

% this method is rigorous only for star shaped surfaces
w = sum(solid_angle(points, faces));

if w<0 && (abs(w)-4*pi)<1000*eps
  s = 'outward';
elseif w>0 && (abs(w)-4*pi)<1000*eps
  s = 'inward';
else
  s = 'unknown';
end

