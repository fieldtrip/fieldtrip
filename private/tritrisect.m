function [l1, l2] = tritrisect(v1, v2, v3, t1, t2, t3)

% TRITRISECT computes the intersection line of a triangle with a plane
% spanned by three vertices v1, v2 and v3.
%
% [l1, l2] = tritrisect(v1, v2, v3, t1, t2, t3)

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: tritrisect.m,v $
% Revision 1.3  2003/03/11 15:35:20  roberto
% converted all files from DOS to UNIX
%
% Revision 1.2  2003/03/04 21:46:20  roberto
% added CVS log entry and synchronized all copyright labels
%

% determine on which side of the plane each vertex lies
p1 = ptriside(v1, v2, v3, t1);
p2 = ptriside(v1, v2, v3, t2);
p3 = ptriside(v1, v2, v3, t3);

if all([p1 p2 p3]==0)
  warning('triangle lies exactly in plane');
  l1 = [nan, nan, nan];
  l2 = [nan, nan, nan];
  return
end

if abs(sum([p1 p2 p3]))==3
  warning('triangle lies on one side of plane');
  l1 = [nan, nan, nan];
  l2 = [nan, nan, nan];
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p1==0 && p2==0 && p3~=0
  % triangle vertex 1 and 2 ly in the plane
  l1 = t1;
  l2 = t2; 
  return
end

if p1==0 && p2~=0 && p3==0
  % triangle vertex 1 and 3 ly in the plane
  l1 = t1;
  l2 = t3; 
  return
end

if p1~=0 && p2==0 && p3==0
  % triangle vertex 2 and 3 ly in the plane
  l1 = t2;
  l2 = t3; 
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p1==0 && p2~=0 && p3~=0
  % triangle vertex 1 lies in the plane
  l1 = t1;
  % triangle edge 2-3 intersects with the plane
  l2 = ltrisect(v1, v2, v3, t2, t3);
  return;
end

if p1~=0 && p2==0 && p3~=0
  % triangle vertex 2 lies in the plane
  l1 = t2;
  % triangle edge 3-1 intersects with the plane
  l2 = ltrisect(v1, v2, v3, t3, t1);
  return;
end

if p1~=0 && p2~=0 && p3==0
  % triangle vertex 3 lies in the plane
  l1 = t3;
  % triangle edge 1-2 intersects with the plane
  l2 = ltrisect(v1, v2, v3, t1, t2);
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p1==p2
  % triangle vertex 3 lies on opposite side
  l1 = ltrisect(v1, v2, v3, t3, t1);
  l2 = ltrisect(v1, v2, v3, t3, t2);
  return;
end

if p2==p3
  % triangle vertex 1 lies on opposite side
  l1 = ltrisect(v1, v2, v3, t1, t2);
  l2 = ltrisect(v1, v2, v3, t1, t3);
  return;
end

if p3==p1
  % triangle vertex 2 lies on opposite side
  l1 = ltrisect(v1, v2, v3, t2, t3);
  l2 = ltrisect(v1, v2, v3, t2, t1);
  return;
end

