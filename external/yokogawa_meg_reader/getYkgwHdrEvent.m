% Get header of the trigger event
% [ Yokogawa MEG Reader toolbox for MATLAB ]
%
% brief:
%   This function retrieves information about trigger event in the specified file.
%
% usage:
%   event = getYkgwHdrEvent(filepath)
%
% arguments:
%   filepath    : file path
%
% return values:
%   event       : structure array of trigger event
%    .sample_no : double, sample number of current event (0 origin)
%    .code      : double, event code (1 origin)
%    .name      : string, event name
%
% rivision history
%   2 : 2011.04.27 : support only EvokedRaw or EvokedAve 
%   1 : 2011.02.14 : 1st argument is modified from file ID to file path.
%   0 : 2010.06.24 : first release
% 
% Copyright (C) 2010-2011 Yokogawa Electric Corporation, All Rights Reserved.