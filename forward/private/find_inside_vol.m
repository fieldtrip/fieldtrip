function [inside, outside] = find_inside_vol(pos, vol);

% FIND_INSIDE_VOL locates dipole locations inside/outside the source
% compartment of a volume conductor model.
% 
% [inside, outside] = find_inside_vol(pos, vol)
%
% This function is obsolete and its use in other functions should be replaced 
% by inside_vol

% Copyright (C) 2003-2007, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information


inside  = inside_vol(pos, vol);
% replace boolean vector with indexing vectors
outside = find(~inside);
inside  = find(inside);
