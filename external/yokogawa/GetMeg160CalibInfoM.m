%  Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
% 
%  usage:
%    calib_info = GetMeg160CalibInfoM( fid )
% 
%  arguments:
%    fid             : file ID
% 
%  return values:
%    calib_info      : (m x 3) matrix of calibration information. (m: channel_count)
%                        detail of n: (If the sensor doesn't have following item, the item's entry will be 'Inf'.)
%                            1   : channel number in Meg160(zero start)
%                            2   : gain [Tesla/V]
%                            3   : offset voltage [V]
%  
%  confirmation of revision:
%   GetMeg160CalibInfoM( Inf ) will show and return revision of this function.
%
