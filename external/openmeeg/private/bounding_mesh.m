function [inside] = bounding_mesh(pos, pnt, tri);

% BOUNDING_MESH determines if a point is inside/outside a triangle mesh 
% whereby the bounding triangle mesh should be closed.
%
% [inside] = bounding_mesh(pos, pnt, tri)
%
% where
%   pos     position of point of interest (can be 1x3 or Nx3)
%   pnt     bounding mesh vertices
%   tri     bounding mesh triangles
%
% See also SOLID_ANGLE

% Copyright (C) 2003, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

global fb;
if isempty(fb)
  fb = 0;
end

npos = size(pos, 1);
npnt = size(pnt, 1);
ntri = size(tri, 1);

% determine a cube that encompases the boundary triangulation
bound_min = min(pnt);
bound_max = max(pnt);

% determine a sphere that is completely inside the boundary triangulation
bound_org = mean(pnt);
bound_rad = sqrt(min(sum((pnt - repmat(bound_org, size(pnt,1), 1)).^2, 2)));

inside = zeros(npos, 1);
for i=1:npos
  if fb
    fprintf('%6.2f%%', 100*i/npos);
  end
  if any(pos(i,:)<bound_min) || any(pos(i,:)>bound_max)
    % the point is outside the bounding cube
    inside(i) = 0;
    if fb, fprintf(' outside the bounding cube\n'); end
  elseif sqrt(sum((pos(i,:)-bound_org).^2, 2))<bound_rad 
    % the point is inside the interior sphere
    inside(i) = 1;
    if fb, fprintf(' inside the interior sphere\n'); end
  else
    % the point is inside the bounding cube but outside the interior sphere
    % compute the total solid angle of the surface, which is zero for a point outside
    % the triangulation and 4*pi or -4*pi for a point inside (depending on the triangle
    % orientation)
    tmp = pnt - repmat(pos(i,:), npnt, 1);
    solang = solid_angle(tmp, tri);
    if any(isnan(solang))
      inside(i) = nan;
    elseif (abs(sum(solang))-2*pi)<0
      % total solid angle is (approximately) zero
      inside(i) = 0;
    elseif (abs(sum(solang))-2*pi)>0
      % total solid angle is (approximately) plus or minus 4*pi
      inside(i) = 1;
    end
    if fb, fprintf(' solid angle\n'); end
  end
end

