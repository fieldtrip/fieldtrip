% Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
%
% usage:
%   [trigger_event] = GetMeg160TriggerEventM( fid )
%
% arguments:
%   fid             : file ID
%
% return values:
%   trigger_event   : n x 1 vector which include multi trigger code (n: epoch count of evoked averaged data)
%
function    [trigger_event] = GetMeg160TriggerEventM( fid )