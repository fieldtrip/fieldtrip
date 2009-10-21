function [volume] = normalisevolume(cfg, volume)

% VOLUMENORMALISE normalizes anatomical and functional data
% to a template anatomical MRI
%
% Warning: This function is deprecated, it has been renamed to VOLUMENORMALISE
% Warning: backward compatibility will be removed in the future

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% $Log: normalisevolume.m,v $
% Revision 1.15  2007/01/09 09:52:07  roboos
% changed the warning a bit
%
% Revision 1.14  2006/03/14 08:09:22  roboos
% added copyrigth and cvs log statement
% 

warning('This function is deprecated, it has been renamed to VOLUMENORMALISE');
warning('backward compatibility will be removed in the future');

[volume] = volumenormalise(cfg, volume);
