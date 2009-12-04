% Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
%
% usage:
%   [system_id version revision system_name] = GetMeg160SystemInfoM( fid )
%
% arguments:
%   fid             : file ID
%
% return values:
%   system_id   : system id
%   version     : version
%   revision    : revision
%   system_name : system name
function    [system_id, version, revision, system_name]  =   GetMeg160SystemInfoM( fid )