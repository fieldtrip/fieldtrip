classdef MultiscaleElectrophysiologyData < handle
    % Class MULTISCALEELECTROPHYSIOLOGYDATA process MED channel data

    % Copyright 2023 Richard J. Cui. Created: Mon 01/30/2023 10:01:06.104 PM
    % $Revision: 0.4 $  $Date: Tue 08/01/2023 12:27:18.039 AM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % =====================================================================
    % properties
    % =====================================================================
    % MED information
    % ---------------
    properties (SetAccess = protected)
        MEDVersion (1, 1) double = NaN % MED version
        MPS = 1e6 % microseconds per seconds
    end % properties

    % MED channel information
    % -----------------------
    properties (SetAccess = protected, Hidden = true)
        FilePath % [str] filepath of MED channel file
        FileName % [str] filename of MED channel file including ext
        Continuity % [table] data segments of conituous sampling (see
        % analyzeContinuity.m for the detail)
        ContinuityCorrected % [table] data segments of conituous sampling corrected
        ChanSamplingFreq % sampling frequency of channel (Hz)
        SampleTimeInterval % sample time interval = [lower, upper] (uUTC),
        % indicating the lower and upper bound of the time interval between
        % two successive samples
    end % properties

    methods

        function this = MultiscaleElectrophysiologyData()

        end

    end % methods

    % other methods
    % -------------
    methods
        varargout = getSampleTimeInterval(this, varargin) % bound of sampling interval
        varargout = SampleTime2Index(this, varargin) % time --> index
        varargout = SampleIndex2Time(this, varargin) % index --> time
    end % methods

end % classdef

% [EOF]
