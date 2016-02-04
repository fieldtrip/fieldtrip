function [proj, dist] = plinprojn(l1,l2,r,flag)

% PLINPROJN projects a point onto a line or linepiece
%
% [proj, dist] = plinprojn(l1, l2, r, flag)
% 
% where l1 and l2 are Nx3 matrices with the begin and endpoints of the linepieces, 
% and r is the point that is projected onto the lines
% This is a vectorized version of Robert's ptriproj function and is
% generally faster than a for-loop around the mex-file.
%
% the optional flag can be:
%   0 (default)  project the point anywhere on the complete line
%   1            project the point within or on the edge of the linepiece

if nargin<4
  flag = false;
end

v  = l2-l1;                  % vector from l1 to l2
dp = bsxfun(@minus, r, l1);  % vector from l1 to r
t  = sum(dp.*v,2)./sum(v.^2,2);

if flag,
  t(t<0) = 0;
  t(t>1) = 1;
end

proj = l1 + [t.*v(:,1) t.*v(:,2) t.*v(:,3)];
dist = sqrt((r(:,1)-proj(:,1)).^2 + (r(:,2)-proj(:,2)).^2 + (r(:,3)-proj(:,3)).^2);
