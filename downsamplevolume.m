function [volume] = downsamplevolume(cfg, volume);

% VOLUMEDOWNSAMPLE downsamples an anatomical MRI or source reconstruction
% and optionally normalizes its coordinate axes, keeping the homogenous
% transformation matrix correct.
%
% Warning: This function is deprecated, it has been renamed to VOLUMEDOWNSAMPLE
% Warning: backward compatibility will be removed in the future

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% $Log: downsamplevolume.m,v $
% Revision 1.12  2007/01/09 09:52:07  roboos
% changed the warning a bit
%
% Revision 1.11  2006/03/14 08:09:22  roboos
% added copyrigth and cvs log statement
% 

warning('This function is deprecated, it has been renamed to VOLUMEDOWNSAMPLE');
warning('backward compatibility will be removed in the future');

[volume] = volumedownsample(cfg, volume);
