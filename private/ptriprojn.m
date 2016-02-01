function [proj, dist] = ptriprojn(v1, v2, v3, r, flag)

% PTRIPROJN projects a point onto the plane going through a set of
% triangles
%
% Use as
%   [proj, dist] = ptriprojn(v1, v2, v3, r, flag)
% where v1, v2 and v3 are Nx3 matrices with vertex positions of the triangles, 
% and r is the point that is projected onto the planes spanned by the vertices
% This is a vectorized version of Robert's ptriproj function and is
% generally faster than a for-loop around the mex-file.
%
% the optional flag can be:
%   0 (default)  project the point anywhere on the complete plane
%   1            project the point within or on the edge of the triangle

if nargin<5
  flag = false;
end

[la, mu, dist, proj] = lmoutrn(v1, v2, v3, r);
switch flag
  case 0
    % finished
  case 1
    % for each of the six regions outside the triangle the projection
    % point is on the edge, or on one of the corners
    
    % la<0 & mu>=0
    sel= la<0;
    [proj(sel,:), dist(sel)] = plinprojn(v1(sel,:), v3(sel,:), r, 1);
    
    % mu<0 &
    sel = mu<0 & la>=0;
    [proj(sel,:), dist(sel)] = plinprojn(v1(sel,:), v2(sel,:), r, 1);
    
    % la+mu>1 & mu>0 & la>0 -> project onto vec2
    sel = la+mu>1 & mu>=0 & la>=0;
    [proj(sel,:), dist(sel)] = plinprojn(v2(sel,:), v3(sel,:), r, 1);
    
end
