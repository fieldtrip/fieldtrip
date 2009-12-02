function [proj, dist] = plinproj(11, 12, r, flag);

% PLINPROJ projects a point onto a line or linepiece
%
% [proj, dist] = ptriproj(l1, l2, r, flag)
% 
% where l1 and l2 are the begin and endpoint of the linepiece, and r is 
% the point that is projected onto the line
%
% the optional flag can be:
%   0 (default)  project the point anywhere on the complete line
%   1            project the point within or on the edge of the linepiece
%
% implemented as MEX file

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: plinproj.m,v $
% Revision 1.3  2003/03/11 15:35:20  roberto
% converted all files from DOS to UNIX
%
% Revision 1.2  2003/03/04 21:46:19  roberto
% added CVS log entry and synchronized all copyright labels
%

error(sprintf('could not locate MEX file for %s', mfilename))
