function [proj, dist] = ptriproj(v1, v2, v3, r, flag);

% PTRIPROJ projects a point onto the plane going through a triangle
%
% [proj, dist] = ptriproj(v1, v2, v3, r, flag)
% 
% where v1, v2 and v3 are three vertices of the triangle, and r is 
% the point that is projected onto the plane spanned by the vertices
%
% the optional flag can be:
%   0 (default)  project the point anywhere on the complete plane
%   1            project the point within or on the edge of the triangle
%
% implemented as MEX file

% Copyright (C) 2002, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

error(sprintf('could not locate MEX file for %s', mfilename))
