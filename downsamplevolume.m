function [volume] = downsamplevolume(cfg, volume);

% VOLUMEDOWNSAMPLE downsamples an anatomical MRI or source reconstruction
% and optionally normalizes its coordinate axes, keeping the homogenous
% transformation matrix correct.
%
% Warning: This function is deprecated, it has been renamed to VOLUMEDOWNSAMPLE
% Warning: backward compatibility will be removed in the future

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

warning('This function is deprecated, it has been renamed to VOLUMEDOWNSAMPLE');
warning('backward compatibility will be removed in the future');

[volume] = volumedownsample(cfg, volume);
