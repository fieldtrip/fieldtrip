%  Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
% 
%  usage:
%    data = GetMeg160ContinuousRawDataM( fid, start_sample, sample_length )
% 
%  arguments:
%    fid             : file ID
%    start_sample    : start sample number for retrieving data. (zero start)
%    sample_length   : sample length for retrieving data.
% 
%  return values:
%    data            : (m x n) matrix of measured data [AD] (class int16)
%                        m: channel_count
%                        n: 1+sample_length (n=1: channel number in Meg160(zero start) )
%  
%  confirmation of revision:
%   GetMeg160ContinuousRawDataM( Inf ) will show and return revision of this function.
%
