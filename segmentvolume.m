function [segment] = segmentvolume(cfg, mri)

% VOLUMESEGMENT segments an anatomical MRI into gray matter, white
% matter, and CSF compartments.
%
% Warning: this function is deprecated, it has been renamed to VOLUMESEGMENT
% Warning: backward compatibility will be removed in the future

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% $Log: segmentvolume.m,v $
% Revision 1.13  2007/01/09 09:52:07  roboos
% changed the warning a bit
%
% Revision 1.12  2006/03/14 08:09:22  roboos
% added copyrigth and cvs log statement
% 

warning('This function is deprecated, it has been renamed to VOLUMESEGMENT');
warning('backward compatibility will be removed in the future');

[volume] = volumesegment(cfg, volume)
