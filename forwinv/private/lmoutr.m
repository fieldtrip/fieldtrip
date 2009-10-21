function [la, mu, dist] = lmoutr(v1, v2, v3, r);

% LMOUTR computes the la/mu parameters of a point projected to a triangle
%
% [la, mu, dist] = lmoutr(v1, v2, v3, r)
% 
% where v1, v2 and v3 are three vertices of the triangle, and r is 
% the point that is projected onto the plane spanned by the vertices
%
% also implemented as MEX file

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: lmoutr.m,v $
% Revision 1.4  2007/01/03 17:02:26  roboos
% removed an incorrect statement from the help
%
% Revision 1.3  2003/03/11 15:35:20  roberto
% converted all files from DOS to UNIX
%
% Revision 1.2  2003/03/04 21:46:19  roberto
% added CVS log entry and synchronized all copyright labels
%

% compute la/mu parameters
vec0 = r  - v1;
vec1 = v2 - v1;
vec2 = v3 - v2;
vec3 = v3 - v1;

tmp   = [vec1' vec3'] \ (vec0');
la    = tmp(1);
mu    = tmp(2);

% determine the projection onto the plane of the triangle
proj  = v1 + la*vec1 + mu*vec3;

% determine the distance from the original point to its projection
dist = norm(r-proj);

