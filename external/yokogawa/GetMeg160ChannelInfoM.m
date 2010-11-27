%  Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
% 
%  usage:
%    channel_info = GetMeg160ChannelInfoM( fid )
% 
%  arguments:
%    fid             : file ID
% 
%  return values:
%    channel_info    : (m x n) matrix of channel information. (m: channel_count, n: 11)
%                        detail of n: (If the sensor doesn't have following item, the item's entry will be 'Inf'.)
%                            1   : channel number in Meg160(zero start)
%                            2   : sensor type ID (ref. table in below)
%                            3-5 : sensor position(x,y,z) [m]
%                            6-9 : sensor direction(zdir1, xdir1, zdir2, xdir2) [deg]
%                            10  : coil size [m]
%                            11  : baseline [m]
% 
%    < table: sensor type ID >
%        NullChannel         = 0;
%        MagnetoMeter        = 1;
%        AxialGradioMeter    = 2;
%        PlannerGradioMeter  = 3;
%        RefferenceChannelMark           = hex2dec('0100');
%        RefferenceMagnetoMeter          = bitor( RefferenceChannelMark, MagnetoMeter );
%        RefferenceAxialGradioMeter      = bitor( RefferenceChannelMark, AxialGradioMeter );
%        RefferencePlannerGradioMeter    = bitor( RefferenceChannelMark, PlannerGradioMeter);
%        TriggerChannel  = -1;
%        EegChannel      = -2;
%        EcgChannel      = -3;
%        EtcChannel      = -4;
%  
%  confirmation of revision:
%   GetMeg160ChannelInfoM( Inf ) will show and return revision of this function.
