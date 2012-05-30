% Get header of the system information
% [ Yokogawa MEG Reader toolbox for MATLAB ]
%
% brief:
%   This function retrieves information of the system in the file specified file.
%
% usage:
%   system_info = getYkgwHdrSystem(filepath)
%
% arguments:
%   filepath                  : file path
%
% return values:
%   system_info : structure of system information
%      .version               : double, system version
%      .revision              : double, system revision
%      .system_id             : double, system ID
%      .system_name           : string, system name
%      .model_name            : string, model name
%
% rivision history
%   1 : 2011.02.14 : 1st argument is modified from file ID to file path.
%   0 : 2010.06.24 : first release
% 
% Copyright (C) 2010-2011 Yokogawa Electric Corporation, All Rights Reserved.
