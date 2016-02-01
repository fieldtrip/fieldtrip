% Get header of the digitization information
% [ Yokogawa MEG Reader toolbox for MATLAB ]
%
% brief:
%   This function retrieves information of the digitization in the specified file.
%
% usage:
%   digitize = getYkgwHdrDigitize(filepath)
%
% arguments:
%   filepath              : file path
%
% return values:
%   digitize
%    .info : structure of system information
%      .digitizer_file    : string, file path of digitizer file
%      .done              :   bool, is matching done?
%      .meg2digitizer     : 4 x 4 matrix, to transform MEG coordinate to Digitizer coordinate 
%      .digitizer2meg     : 4 x 4 matrix, to transform Digitizer coordinate to MEG coordinate 
%    .point               : structure array of point data
%      .name              : string, point name
%      .x                 : double, x-coordinate on digitizer coordinate [meter]
%      .y                 : double, y-coordinate on digitizer coordinate [meter]
%      .z                 : double, z-coordinate on digitizer coordinate [meter]
%
% rivision history
%   2 : 2011.04.27 : fix .meg2digitizer .digitizer2meg
%   1 : 2011.02.14 : 1st argument is modified from file ID to file path.
%   0 : 2010.06.24 : first release
% 
% Copyright (C) 2010-2011 Yokogawa Electric Corporation, All Rights Reserved.
