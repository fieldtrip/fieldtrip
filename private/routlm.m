function [proj] = routlm(v1, v2, v3, la, mu);

% ROUTLM computes the projection of a point from its la/mu parameters
% these equal the "Barycentric" coordinates
%
% [proj] = routlm(v1, v2, v3, la, mu)
% 
% where v1, v2 and v3 are three vertices of the triangle
%
% also implemented as MEX file

% Copyright (C) 2002, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% determine the projection onto the plane of the triangle
proj  = (1-la-mu)*v1 + la*v2 + mu*v3;
