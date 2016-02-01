function [side] = ptriside(v1, v2, v3, r, tolerance)

% PTRISIDE determines the side of a plane on which a set of points lie. it
% returns 0 for the points that lie on the plane
%
% [side] = ptriside(v1, v2, v3, r)
% 
% the side of points r is determined relative to the plane spanned by
% vertices v1, v2 and v3. v1,v2 and v3 should be 1x3 vectors. r should be a
% Nx3 matrix

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: ptriside.m,v $
% Revision 1.3  2003/03/11 15:35:20  roberto
% converted all files from DOS to UNIX
%
% Revision 1.2  2003/03/04 21:46:19  roberto
% added CVS log entry and synchronized all copyright labels
%

if nargin<5
  tolerance = 100*eps;
end

n = size(r,1);
a = r  - ones(n,1)*v1;
b = v2 - v1;
c = v3 - v1;
d = crossproduct(b, c);
val = dotproduct(a, d);

side = zeros(n, 1);
side(val >  tolerance) =  1;
side(val < -tolerance) = -1;

%if val>tolerance
%  side=1;
%elseif val<-tolerance
%  side=-1;
%else
%  side=0;
%end

% subfunction without overhead to speed up
function c = crossproduct(a, b)

c(1) = a(2)*b(3)-a(3)*b(2);
c(2) = a(3)*b(1)-a(1)*b(3);
c(3) = a(1)*b(2)-a(2)*b(1);

% subfunction without overhead to speed up, input a can be a matrix
function d = dotproduct(a, b)

d = a(:,1)*b(1)+a(:,2)*b(2)+a(:,3)*b(3);
