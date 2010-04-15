function ft_sourcewrite(cfg, volume)

% FT_VOLUMEWRITE exports source analysis results to an Analyze MRI file
% that can subsequently be read into BrainVoyager or MRIcro
%
% Warning: FT_NORMALISEVOLUME has been renamed to FT_VOLUMENORMALISE
% Warning: backward compatibility will be removed in the future

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

warning('SOURCEWRITE has been renamed to VOLUMEWRITE');
warning('backward compatibility will be removed in the future');

ft_volumewrite(cfg, volume);
