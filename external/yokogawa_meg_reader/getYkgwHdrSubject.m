% Get header of the subject information
% [ Yokogawa MEG Reader toolbox for MATLAB ]
%
% brief:
%   This function retrieves information of the subject in the specified file.
%
% usage:
%   subject = getYkgwHdrSubject(filepath)
%
% arguments:
%   filepath     : file path
%
% return values:
%   subject      : structure of subject information
%     .id        : string, id
%     .name      : string, name
%     .birthday  : string, birthday
%     .sex       : string, sex
%     .handed    : string, handed
%
% rivision history
%   1 : 2011.02.14 : 1st argument is modified from file ID to file path.
%   0 : 2010.06.24 : first release
% 
% Copyright (C) 2010-2011 Yokogawa Electric Corporation, All Rights Reserved.
