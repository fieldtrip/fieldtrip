function [segment] = segmentvolume(cfg, mri)

% VOLUMESEGMENT segments an anatomical MRI into gray matter, white
% matter, and CSF compartments.
%
% Warning: this function is deprecated, it has been renamed to VOLUMESEGMENT
% Warning: backward compatibility will be removed in the future

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

warning('This function is deprecated, it has been renamed to VOLUMESEGMENT');
warning('backward compatibility will be removed in the future');

[volume] = volumesegment(cfg, volume)
