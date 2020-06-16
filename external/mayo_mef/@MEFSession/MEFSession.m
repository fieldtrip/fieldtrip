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
    % $Revision: 0.5 $  $Date: Thu 04/02/2020  2:13:08.181 PM $
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
        SelectedChannel     % channels selected
        StartEnd            % start and end points to import the session
        SEUnit              % unit of StartEnd
    end % properties

    % properties of session information
    % ---------------------------------
    properties        
        SessionPath         % session directory
        Password            % password structure of the session
        ChannelName         % channel names
        SamplingFrequency   % in Hz
        Samples             % number of samples
        DataBlocks          % number of data blocks
        TimeGaps            % number of discountinuity time gaps
        BeginStop           % Begin and stop indexes of entire signal
        Unit                % unit of BeginStop
        Institution         % name of the institute
        SubjectID           % identification of the subject
        AcquisitionSystem   % name of the system to record the session
        CompressionAlgorithm % name of compression algorithm
        SessionInformation  % table of session information (see get_sessinfo.m)
        SessionContinuity   % [table] data segments of conituous sampling (see
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
        varargout = get_sessinfo(this) % get sess info from data
        [X, t] = importSession(this, varargin) % import a session
        record_offset = getSessionRecordOffset(this, varargin) % get offset time of recording in specified unit
        rel_time = abs2relativeTimePoint(this, abs_time, unit) % absolute to relative time points
        abs_time = relative2absTimePoint(this, rel_time, unit) % relative to absolute time points
        out_time = SessionUnitConvert(this, in_time, varargin) % convert units of relative time points
        [path_to_sess, sess_name, sess_ext] = get_sess_parts(this, varargin) % get the parts of session path
    end % methods
end

% [EOF]