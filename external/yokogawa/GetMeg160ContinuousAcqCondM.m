% Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
%
% usage:
%   [sample_rate, sample_count] = GetMeg160ContinuousAcqCondM( fid )
%
% arguments:
%   fid             : file ID
%
% return values:
%   sample_rate     : sample rate [Hz]
%   sample_count    : actual sampled count [sample]
function    [sample_rate, sample_count]  =   GetMeg160ContinuousAcqCondM( fid )