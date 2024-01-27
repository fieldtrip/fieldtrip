classdef MEFSession < handle
    % Class MEFSESSION processes MEF session
    %
    % Syntax:
    %   this = MEFSession;
    %
    % Input(s):
    %
    % Output(s):
    %
    % See also .

    % Copyright 2020 Richard J. Cui. Created: Thu 02/06/2020 10:07:26.965 AM
    % $Revision: 0.6 $  $Date: Wed 09/13/2023 10:48:18.575 PM $
    %
    % 1026 Rocky Creek Dr NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % =====================================================================
    % properties
    % =====================================================================
    % properties of importing session
    % -------------------------------
    properties
        SelectedChannel % channels selected
        StartEnd % start and end points to import the session
        SEUnit % unit of StartEnd
    end % properties

    % properties of session information
    % ---------------------------------
    properties
        SessionPath % session directory
        Password % password structure of the session
        ChannelName % channel names
        SamplingFrequency % in Hz
        Samples % number of samples
        DataBlocks % number of data blocks
        TimeGaps % number of discountinuity time gaps
        BeginStop % Begin and stop indexes of entire signal
        Unit % unit of BeginStop
        Institution % name of the institute
        SubjectID % identification of the subject
        AcquisitionSystem % name of the system to record the session
        CompressionAlgorithm % name of compression algorithm
        SessionInformation % table of session information (see get_sessinfo.m)
        SessionContinuity % [table] data segments of conituous sampling (see
        % analyzeContinuity.m for the detail)
    end % properties

    % =====================================================================
    % methods
    % =====================================================================
    % the constructor
    % ----------------
    methods

        function this = MEFSession()

        end % function

    end % methods

    % other metheds
    % -------------
    methods
        varargout = abs2relativeTimePoint(this, varargin) % absolute to relative time points
        varargout = get_sessinfo(this, varargin) % get sess info from data
        varargout = get_sess_parts(this, varargin) % get the parts of session path
        varargout = getSessionRecordOffset(this, varargin) % get offset time of recording in specified unit
        varargout = importSession(this, varargin) % import a session
        varargout = relative2absTimePoint(this, varargin) % relative to absolute time points
        varargout = SessionUnitConvert(this, varargin) % convert units of relative time points
    end % methods

end

% [EOF]
