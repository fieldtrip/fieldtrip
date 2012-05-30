% Get measured data
% [ Yokogawa MEG Reader toolbox for MATLAB ]
%
% brief:
%   This function retrieves the measurement data of whole channel by the specified file and sample range.
%
% usage:
%   data = getYkgwData(filepath, start_sample, sample_length)
%
% arguments:
%   filepath        : file path
%
%   start_sample    : The start number corresponding to each acquisition type is as follows :
%    Continuous Raw : Start sample number for retrieving data. (0 origin)
%    Evoked Average : Start sample number for retrieving data. (0 origin)
%    Evoked Raw     : Start frame  number for retrieving data. (1 origin)
%      note : When both start_sample and sample_length are omitted, you can get data of whole samples.
%
%   sample_length   : The number of samples or trials corresponding to each acquisition type is as follows :
%    Continuous Raw : Number of samples for retrieving data.
%    Evoked Average : Number of frames for retrieving data.
%    Evoked Raw     : Number of frames for retrieving data.
%      note : When this parameter is omitted or is specified as 'Inf',
%              you can get data from start_sample to the end of sample(frame).
%
% return values:
%   data            : (number of channel x number of sample) double matrix of measured data 
%          The unit corresponding to each channel type is as follows :
%           AxialGradioMeter              [Tesla]
%           PlannerGradioMeter            [Tesla]
%           MagnetoMeter                  [Tesla]
%           RefferenceMagnetoMeter        [Tesla]
%           RefferenceAxialGradioMeter    [Tesla]
%           RefferencePlannerGradioMeter  [Tesla]
%           TriggerChannel                [Volt]
%           EegChannel                    [Volt] *This has already been reflected EEG gain
%           EcgChannel                    [Volt] *This has already been reflected ECG gain
%           EtcChannel                    [Volt]
%           NullChannel                   [Volt]
%
%
% rivision history
%   4 : 2011.03.24 : Support GetSqf changes
%   3 : 2011.03.03 : Support MegLaboratory R1.4.8 calibration for non-MEG
%                    channels whose sensitivities are not specified.
%   2 : 2011.02.14 : Fix EvokedRaw error if sample_length argument is omitted.
%                    1st argument is modified from file ID to file path.
%   1 : 2011.01.24 : Not apply averaged gain to reference channels if system id is 1000 later. 
%   0 : 2010.06.24 : first release
% 
% Copyright (C) 2010-2011 Yokogawa Electric Corporation, All Rights Reserved.
