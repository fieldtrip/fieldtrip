function [X, Y, Z] = intersect_plane(pnt, dhk, v1, v2, v3)

% INTERSECT_PLANE intersection between a triangulated surface and a plane
% it returns the coordinates of the vertices which form a contour
%
% [X, Y, Z] = intersect_plane(pnt, dhk, v1, v2, v3)
%
% the intersecting plane is spanned by the vertices v1, v2, v3

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: intersect_plane.m,v $
% Revision 1.3  2003/03/11 15:35:20  roberto
% converted all files from DOS to UNIX
%
% Revision 1.2  2003/03/04 21:46:19  roberto
% added CVS log entry and synchronized all copyright labels
%

npnt = size(pnt,1);
ndhk = size(dhk,1);
side = zeros(npnt,1);
for i=1:npnt
  side(i) = ptriside(v1, v2, v3, pnt(i,:));
end

% find the triangles which have vertices on both sides of the plane
indx = find(abs(sum(side(dhk),2))~=3);
cnt1 = zeros(length(indx), 3);
cnt2 = zeros(length(indx), 3);

for i=1:length(indx)
  cur = dhk(indx(i),:);
  tmp = side(cur);
  l1 = pnt(cur(1),:);
  l2 = pnt(cur(2),:);
  l3 = pnt(cur(3),:);
  if tmp(1)==tmp(2)
    cnt1(i,:) = ltrisect(v1, v2, v3, l3, l1);
    cnt2(i,:) = ltrisect(v1, v2, v3, l3, l2);
  elseif tmp(1)==tmp(3)
    cnt1(i,:) = ltrisect(v1, v2, v3, l2, l1);
    cnt2(i,:) = ltrisect(v1, v2, v3, l2, l3);
  elseif tmp(2)==tmp(3)
    cnt1(i,:) = ltrisect(v1, v2, v3, l1, l2);
    cnt2(i,:) = ltrisect(v1, v2, v3, l1, l3);
  elseif tmp(1)==0 && tmp(2)==0
    cnt1(i,:) = l1;
    cnt2(i,:) = l2;
  elseif tmp(1)==0 && tmp(3)==0
    cnt1(i,:) = l1;
    cnt2(i,:) = l3;
  elseif tmp(2)==0 && tmp(3)==0
    cnt1(i,:) = l2;
    cnt2(i,:) = l3;
  elseif tmp(1)==0 && tmp(2)~=tmp(3)
    cnt1(i,:) = l1;
    cnt2(i,:) = ltrisect(v1, v2, v3, l2, l3);
  elseif tmp(2)==0 && tmp(3)~=tmp(1)
    cnt1(i,:) = l2;
    cnt2(i,:) = ltrisect(v1, v2, v3, l3, l1);
  elseif tmp(3)==0 && tmp(1)~=tmp(2)
    cnt1(i,:) = l3;
    cnt2(i,:) = ltrisect(v1, v2, v3, l1, l2);
  elseif tmp(1)==0
    cnt1(i,:) = l1;
    cnt2(i,:) = l1;
  elseif tmp(2)==0
    cnt1(i,:) = l2;
    cnt2(i,:) = l2;
  elseif tmp(3)==0
    cnt1(i,:) = l3;
    cnt2(i,:) = l3;
  end
end

X = [cnt1(:,1) cnt2(:,1)];
Y = [cnt1(:,2) cnt2(:,2)];
Z = [cnt1(:,3) cnt2(:,3)];

