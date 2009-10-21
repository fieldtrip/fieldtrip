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
% $Log: routlm.m,v $
% Revision 1.4  2003/03/11 15:35:20  roberto
% converted all files from DOS to UNIX
%
% Revision 1.3  2003/03/04 21:46:19  roberto
% added CVS log entry and synchronized all copyright labels
%

% determine the projection onto the plane of the triangle
proj  = (1-la-mu)*v1 + la*v2 + mu*v3;
