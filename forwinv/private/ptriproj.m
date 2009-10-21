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
% $Log: ptriproj.m,v $
% Revision 1.4  2003/03/19 12:22:23  roberto
% fixed bug in function definition (added flag input parameter)
%
% Revision 1.3  2003/03/11 15:35:20  roberto
% converted all files from DOS to UNIX
%
% Revision 1.2  2003/03/04 21:46:19  roberto
% added CVS log entry and synchronized all copyright labels
%

error(sprintf('could not locate MEX file for %s', mfilename))
