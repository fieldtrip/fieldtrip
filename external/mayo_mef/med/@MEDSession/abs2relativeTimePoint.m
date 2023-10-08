function rel_time = abs2relativeTimePoint(this, abs_time, unit)
    % MEDSESSION.ABS2RELATIVETIMEPOINT convert absolute time point to relative one
    %
    % Syntax:
    %   rel_time = abs2relativeTimePoint(this, abs_time, unit)
    %
    % Input(s):
    %   this            - [obj] MEDSession object
    %   abs_time        - [array] absolute time points
    %   unit            - [char] unit of the time points
    %
    % Output(s):
    %   rel_time        - [array] relative time points, which are relative to
    %                     the beginning of recording
    %
    % Note:
    %
    % Seaa slso .

    % Copyright 2023 Richard J. Cui. Created: Wed 09/13/2023 10:48:18.575 PM
    % $Revision: 0.1 $  $Date: Wed 09/13/2023 10:48:18.575 PM $
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
        abs_time (:, 1) double
        unit (1, 1) string ...
            {mustBeMember(unit, ["index", "uutc", "msec", "second", "minute", "hour", "day"])}
    end % positional

    % =========================================================================
    % main
    % =========================================================================
    offset = this.getSessionRecordOffset(unit);
    rel_time = abs_time - offset;
end

% =========================================================================
% subroutines
% =========================================================================

% [EOF]
