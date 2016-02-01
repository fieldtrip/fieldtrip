% Get header of the system information
% [ Yokogawa MEG Reader toolbox for MATLAB ]
%
% brief:
%   This function retrieves information about data acquisition condition in the specified file.
%
% usage:
%   acq_cond = getYkgwHdrAcqCond(filepath)
%
% arguments:
%   filepath                        : file path
%
% return values:
%   acq_cond : structure of acquisition condition
%    .acq_type                      : double, acquisition type (1 : ContinuousRaw, 2 : EvokedAverage, 3 : EvokedRaw)
%         AcqTypeContinuousRaw      = 1;
%         AcqTypeEvokedAve          = 2;
%         AcqTypeEvokedRaw          = 3;
%
%   [when ContinuousRaw(*.con) acquisition type]
%    .sample_rate                   : double, sampling rate [Hz]
%    .sample_count                  : double, number of sample which actually acquired [sample]
%    .specified_sample_count        : double, The number of samples which were specified before starting acquisition [sample]
%
%   [when EvokedRaw(*.raw) or EvokedAve(*.ave) acquisition type] 
%    .sample_rate                   : double, sampling rate [Hz]
%    .frame_length                  : double, frame length [sample]
%    .pretrigger_length             : double, pretrigger length [sample]
%    .average_count                 : double, The number of trials(frames) which were actually acquired [trial]
%    .specified_average_count       : double, The number of trials(frames) which were specified before starting acquisition [trials]
%
%    .multi_trigger                 : The structure of multi trigger information.
%       .enable                     : boolean, Is multi trigger mode ? (true : multi trigger mode)
%       .count                      : double,  Number of multi triggers
%       .list                       : structure array List of multi triggers  
%                                    (If not multi trigger mode, this structure array is set to empty.)
%          .enable                  : boolean, Is current multi trigger set to enable ? (true : enable)
%          .code                    : double,  Event code (1 origin)
%          .name                    : string,  Event name
%          .average_count           : double,  The number of trials(frames) which were actually acquired [trial]
%          .specified_average_count : double,  The number of trials(frames) which were specified before starting acquisition [trials]
%         
% rivision history
%   2 : 2011.04.27 : 'multi_trigger' field was added to output structure.
%   1 : 2011.02.14 : 1st argument is modified from file ID to file path.
%   0 : 2010.06.24 : first release
% 
% Copyright (C) 2010-2011 Yokogawa Electric Corporation, All Rights Reserved.
