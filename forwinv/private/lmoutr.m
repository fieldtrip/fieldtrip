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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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

