function [volume] = normalisevolume(cfg, volume)

% VOLUMENORMALISE normalizes anatomical and functional data
% to a template anatomical MRI
%
% Warning: This function is deprecated, it has been renamed to VOLUMENORMALISE
% Warning: backward compatibility will be removed in the future

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

warning('This function is deprecated, it has been renamed to VOLUMENORMALISE');
warning('backward compatibility will be removed in the future');

[volume] = volumenormalise(cfg, volume);
