% Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
%
% usage:
%   acq_type = GetMeg160DataAcqTypeM( fid )
%
% arguments:
%   fid             : file ID
%
% return values:
%   acq_type        : id number of data acquisition type
%                       1: continuous raw data
%                       2: evoked averaged data
%                       3: evoked raw data
function    acq_type  =   GetMeg160DataAcqTypeM( fid )