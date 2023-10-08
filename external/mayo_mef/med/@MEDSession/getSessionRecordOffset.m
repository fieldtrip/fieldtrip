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
    % $Revision: 0.2 $  $Date: Fri 10/06/2023 12:19:28.864 AM $
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
    offset_uutc = this.MetaData.recording_time_offset; % in uutc

    switch unit
        case "index"
            record_offset = 0;
        case "uutc"
            record_offset = offset_uutc;
        case "msec"
            record_offset = offset_uutc / 1e3;
        case "second"
            record_offset = offset_uutc / 1e6;
        case "minute"
            record_offset = offset_uutc / 1e6/60;
        case "hour"
            record_offset = offset_uutc / 1e6/60/60;
        case "day"
            record_offset = offset_uutc / 1e6/60/60/24;
    end % switch

end % funciton

% =========================================================================
% subroutines
% =========================================================================

% [EOF]
