function record_offset = getSessionRecordOffset(this, unit)
    % MEDSESSION.GETSESSIONRECORDOFFSET get offset time of recording in specified unit
    %
    % Syntax:
    %   record_offset = getSessionRecordOffset(this, unit)
    %
    % Input(s):
    %   this            - [obj] MEFSession_2p1 object
    %   unit            - [str] (opt) unit of the output offset (default =
    %                     uUTC)
    %
    % Output(s):
    %   record_offset   - [num] recording offset of time
    %
    % Note:
    %
    % Seaa slso .

    % Copyright 2020 Richard J. Cui. Created: Wed 09/13/2023 11:45:31.702 PM
    % $Revision: 0.3 $  $Date: Tue 10/10/2023 09:54:21.098 PM $
    %
    % 1026 Rocky Creek Dr NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % =========================================================================
    % parse inputs
    % =========================================================================
    arguments
        this (1, 1) MEDSession
        unit (1, 1) string ...
            {mustBeMember(unit, ["index", "uutc", "msec", "second", "minute", "hour", "day"])}
    end % positional

    % =========================================================================
    % main
    % =========================================================================
    start_time = double(this.MetaData.start_time); % in uutc
    start_sample = double(this.MetaData.absolute_start_sample_number); % in sample
    offset_uutc = double(this.MetaData.recording_time_offset); % in uutc

    switch unit
        case "index"
            record_offset = start_sample - 1;
        case "uutc"
            record_offset = offset_uutc + start_time;
        case "msec"
            record_offset = (offset_uutc + start_time) / 1e3;
        case "second"
            record_offset = (offset_uutc +start_time) / 1e6;
        case "minute"
            record_offset = (offset_uutc +start_time) / 1e6/60;
        case "hour"
            record_offset = (offset_uutc +start_time) / 1e6/60/60;
        case "day"
            record_offset = (offset_uutc +start_time) / 1e6/60/60/24;
    end % switch

end % funciton

% =========================================================================
% subroutines
% =========================================================================

% [EOF]
