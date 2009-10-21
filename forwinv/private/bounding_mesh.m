function [inside] = bounding_mesh(pos, pnt, tri);

% BOUNDING_MESH determines if a point is inside/outside a triangle mesh 
% whereby the bounding triangle mesh should be closed.
%
% [inside] = bounding_mesh(pos, pnt, tri)
%
% where
%   pos		position of point of interest (can be 1x3 or Nx3)
%   pnt		bounding mesh vertices
%   tri		bounding mesh triangles
%
% See also SOLID_ANGLE

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: bounding_mesh.m,v $
% Revision 1.11  2006/03/06 09:44:34  roboos
% changed a | into a ||
%
% Revision 1.10  2003/07/24 09:18:17  roberto
% added abs to solid angle to make it invariant for triange orientation
%
% Revision 1.9  2003/03/21 13:40:29  roberto
% fixed small bug (unmatched end)
%
% Revision 1.8  2003/03/21 13:38:46  roberto
% solid angle computation is now available in very fast mex file
% removed projection along x/y/z axes
%
% Revision 1.7  2003/03/21 12:32:35  roberto
% added projection along y and x, similar to z
%
% Revision 1.6  2003/03/21 12:15:23  roberto
% added new approach to determine inside/outside (proj along z)
%
% Revision 1.5  2003/03/21 08:26:46  roberto
% changed to percent (100x)
%
% Revision 1.4  2003/03/21 08:26:08  roberto
% fixed typo
%
% Revision 1.3  2003/03/21 08:25:34  roberto
% small modification to feedback/debugging info
%
% Revision 1.2  2003/03/04 21:35:26  roberto
% added CVS Log keyword
%

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

