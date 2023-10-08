function abs_time = relative2absTimePoint(this, rel_time, unit)
    % MEDSESSION.RELATIVE2ABSTIMEPOINT convert relative time point to absolute one
    %
    % Syntax:
    %   abs_time = abs2relativeTimePoint(this, rel_time, unit)
    %
    % Input(s):
    %   this            - [obj] MEFSession_2p1 object
    %   rel_time        - [array] relative time points, which are relative to
    %                     the beginning of recording
    %   unit            - [char] unit of the time points
    %
    % Output(s):
    %   abs_time        - [array] absolute time points
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    % See also .

    % Copyright 2023 Richard J. Cui. Created: Sun 09/17/2023 10:40:37.047 PM
    % $Revision: 0.1 $  $Date: Sun 09/17/2023 10:40:37.096 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        this (1, 1) MEDSession
        rel_time (1, :) double
        unit (1, :) char {mustBeMember(unit, {'index', 'uutc', 'msec', 'second', 'minute', 'hour', 'day'})} = 's'
    end % positional

    % ======================================================================
    % main
    % ======================================================================
    offset = this.getSessionRecordOffset(unit);
    abs_time = rel_time + offset;

end % function relative2absTimePoint

% [EOF]
