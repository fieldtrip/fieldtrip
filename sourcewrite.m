function sourcewrite(cfg, volume)

% VOLUMEWRITE exports source analysis results to an Analyze MRI file
% that can subsequently be read into BrainVoyager or MRIcro
%
% Warning: NORMALISEVOLUME has been renamed to VOLUMENORMALISE
% Warning: backward compatibility will be removed in the future

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% $Log: sourcewrite.m,v $
% Revision 1.10  2006/03/14 08:09:22  roboos
% added copyrigth and cvs log statement
% 

warning('SOURCEWRITE has been renamed to VOLUMEWRITE');
warning('backward compatibility will be removed in the future');

volumewrite(cfg, volume);
