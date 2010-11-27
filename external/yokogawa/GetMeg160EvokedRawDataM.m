%  Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
% 
%  usage:
%    data = GetMeg160EvokedRawDataM( fid, start_epoch, epoch_count )
% 
%  arguments:
%    fid             : file ID
%    start_epoch     : start epoch number for retrieving data. (zero start)
%    epoch_count     : epoch count for retrieving data.
% 
%  return values:
%    data            : (m x n) matrix of measured data [AD] (class: int16)
%                        m: channel_count
%                        n: 1+sample_length (n=1: channel number in Meg160(zero start) )
%  
%  confirmation of revision:
%   GetMeg160EvokedRawDataM( Inf ) will show and return revision of this function.
%
