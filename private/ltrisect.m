function [sect] = ltrisect(v1, v2, v3, l1, l2);

% LTRISECT intersects a line with a plane spanned by three vertices
%
% [sect] = ltrisect(v1, v2, v3, l1, l2)
% 
% where v1, v2 and v3 are three vertices spanning the plane, and l1 and l2
% are two points on the line
%
% implemented as MEX file

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: ltrisect.m,v $
% Revision 1.3  2003/03/11 15:35:20  roberto
% converted all files from DOS to UNIX
%
% Revision 1.2  2003/03/04 21:46:19  roberto
% added CVS log entry and synchronized all copyright labels
%

error(sprintf('could not locate MEX file for %s', mfilename))
