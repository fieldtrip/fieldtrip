%  Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
% 
%  usage:
%    [sample_rate, frame_length, pretrigger_length, averaged_count, (multi_trigger_count), (multi_trigger_list)] = GetMeg160EvokedAcqCondM( fid )
% 
%  arguments:
%    fid             : file ID
% 
%  return values:
%    sample_rate         : sample rate [Hz]
%    frame_length        : frame length [sample]
%    pretrigger_length   : pretrigger length [sample]
%    averaged_count      : actual averaged count
%    multi_trigger_count : multi trigger count
%    multi_trigger_list  : Structures of each multi trigger information.
%                            .enable         : flag of enable or diable (enable: true, diable: false)
%                            .attrib         : (reserved for future)
%                            .status         : (reserved for future)
%                            .event_code     : code number of this trigger
%                            .name           : trigger name
%                            .averaged_count : average count
%  
%  confirmation of revision:
%   GetMeg160EvokedAcqCondM( Inf ) will show and return revision of this function.
%
