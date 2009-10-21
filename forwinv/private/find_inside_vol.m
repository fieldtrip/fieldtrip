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
% $Log: find_inside_vol.m,v $
% Revision 1.1  2009/01/21 10:32:38  roboos
% moved from forwinv/* and forwinv/mex/* directory to forwinv/private/* to make the CVS layout consistent with the release version
%
% Revision 1.8  2008/09/20 13:41:35  roboos
% moved content of find_inside_vol to new inside_vol function with slightly different interface
% added wrapper for spm
%


inside  = inside_vol(pos, vol);
% replace boolean vector with indexing vectors
outside = find(~inside);
inside  = find(inside);
