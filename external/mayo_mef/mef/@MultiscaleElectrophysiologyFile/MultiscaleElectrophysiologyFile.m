classdef MultiscaleElectrophysiologyFile < handle
    % Class MULTISCALEELECTROPHYSIOLOGYFILE process MEF channel data
    % 
    % Syntax:
    %   this = MultiscaleElectrophysiologyFile;
    %
    % Input(s):
    %
    % Output(s):
    %
    % Note:
    % 
    % See also .
    
    % Copyright 2020 Richard J. Cui. Created: Tue 02/04/2020  2:21:31.965 PM
    % $Revision: 0.4 $  $Date: Wed 03/11/2020 10:07:12.332 PM $
    %
    % 1026 Rocky Creek Dr NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % =====================================================================
    % properties
    % =====================================================================
    % MEF information
    % ---------------
    properties (SetAccess = protected)
        MEFVersion = []     % MEF version to serve, can be set only in 
                            % constructor
        MPS = 1e6           % microseconds per seconds
    end %  properties: protected
    
    % MEF channel info
    % ----------------
    properties (SetAccess = protected, Hidden = true)
        FilePath            % [str] filepath of MEF channel file
        FileName            % [str] filename of MEF channel file including ext
        Header              % [struct] header information of MEF file
        BlockIndexData      % [table] data of block indices (see 
                            % readBlockIndexData.m for the detail)
        Continuity          % [table] data segments of conituous sampling (see
                            % analyzeContinuity.m for the detail)
        ChanSamplingFreq    % sampling frequency of channel (Hz)
        SampleTimeInterval  % sample time interval = [lower, upper] (uUTC),
                            % indicating the lower and upper bound of the
                            % time interval between two successive samples
    end % properties: protected, hidden
    
    properties
    end % properties

    % =====================================================================
    % methods
    % =====================================================================
    % the constructor
    % ----------------    
    methods
        function this = MultiscaleElectrophysiologyFile()
        
        end % function
    end
    
    % other metheds
    % -------------
    methods
        sti = getSampleTimeInterval(this, varargin) % bound of sampling interval
        [sample_index, sample_yn] = SampleTime2Index(this, varargin) % time --> index
        [sample_time, sample_yn] = SampleIndex2Time(this, varargin) % index --> time
        this = setContinuity(this, cont_table) % set Continuity table
        out_time = SampleUnitConvert(this, in_time, varargin) % convert units of time points
        record_offset = getRecordOffset(this, unit) % get offset time of recording in specified unit
    end % methods
end

% [EOF]