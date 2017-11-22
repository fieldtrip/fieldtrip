function dist = pntdist(p1, p2)

% PNTDIST returns the euclidian distance between two points
%
%  [dist] = pntdist(pnt1, pnt2)
%
% where pnt1 and pnt2 must be Npnt x 3
% or either one can be Npnt x 1

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: pntdist.m,v $
% Revision 1.2  2003/03/04 21:46:19  roberto
% added CVS log entry and synchronized all copyright labels
%

if size(p1,1)==1
  p1 = repmat(p1, size(p2,1), 1);
elseif size(p2,1)==1
  p2 = repmat(p2, size(p1,1), 1);
end

dist = sqrt(sum((p1-p2).^2, 2));

