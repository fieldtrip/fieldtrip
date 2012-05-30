% Get header of the system information
% [ Yokogawa MEG Reader toolbox for MATLAB ]
%
% brief:
%   This function retrieves information about channel in the specified file.
%
% usage:
%   channel_info = getYkgwHdrChannel(filepath)
%
% arguments:
%   filepath           : file path
%
% return values:
%   channel_info       : structure of channel information
%     .channel_count   : double, total number of channel
%     .channel         : structure array ('index 1' corresponds to 'channel 0')
%        .type         : double, channel type as follows:
%          NullChannel                     = 0;
%          MagnetoMeter                    = 1;
%          AxialGradioMeter                = 2;
%          PlanarGradioMeter               = 3;
%          ReferenceMagnetoMeter           = 257;
%          ReferenceAxialGradioMeter       = 258;
%          ReferencePlanarGradioMeter      = 259;
%          TriggerChannel                  = -1;
%          EegChannel                      = -2;
%          EcgChannel                      = -3;
%          EtcChannel                      = -4;
%       .data          : structure of channel specifications, 
%                        fields of each channel type is as follows:
%         [ AxialGradioMeter, ReferenceAxialGradioMeter ]
%           .x        : double, x coordinate [m] of pickup sensor position
%           .y        : double, y coordinate [m] of pickup sensor position
%           .z        : double, z coordinate [m] of pickup sensor position
%           .zdir     : double, sensor direction [degree]
%           .xdir     : double, sensor direction [degree]
%           .baseline : double, length of baseline [m]
%           .size     : double, pickup coil size [m] (diameter of coil circle)
%           .name     : string, abbreviation name
%         [ PlanarGradioMeter, ReferencePlanarGradioMeter ]
%           .x        : double, x coordinate [m] of pickup sensor position
%           .y        : double, y coordinate [m] of pickup sensor position
%           .z        : double, z coordinate [m] of pickup sensor position
%           .zdir1    : double, sensor direction [degree]
%           .xdir1    : double, sensor direction [degree]
%           .zdir2    : double, baseline direction (pickup coil to reference coil) [degree]
%           .xdir2    : double, baseline direction (pickup coil to reference coil) [degree]
%           .baseline : double, length of baseline [m]
%           .size     : double, pickup coil size [m] (a side of square coil)
%         [ MagnetoMeter, ReferenceMagnetoMeter ]
%           .x        : double, x coordinate [m] of coil position
%           .y        : double, y coordinate [m] of coil position
%           .z        : double, z coordinate [m] of coil position
%           .zdir     : double, coil direction [degree]
%           .xdir     : double, coil direction [degree]
%           .size     : double, coil size [m] (diameter of coil circle)
%           .name     : string, abbreviation name
%         [ TriggerChannel, EtcChannel ]
%           .type     : double, type
%           .id       : double, ID
%           .name     : string, abbreviation name
%         [ EegChannel, EcgChannel ]
%           .type     : double, EEG/ECG type
%           .id       : double, EEG/ECG ID
%           .name     : string, abbreviation name
%           .gain     : double, gain
%         [ NullChannel ] no field
%         * All coordinate is based on MEG device coordinate system.
%
% rivision history
%   2 : 2011.04.27 : remove the end of unnecessary spaces for channel name
%   1 : 2011.02.14 : 1st argument is modified from file ID to file path.
%   0 : 2010.06.30 : first release
% 
% Copyright (C) 2010-2011 Yokogawa Electric Corporation, All Rights Reserved.
