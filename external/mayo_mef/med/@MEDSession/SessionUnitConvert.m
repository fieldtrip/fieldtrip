function out_time = SessionUnitConvert(this, in_time, options)
    % MEDSESSION.SESSIONUNITCONVERT convert units of relative time points
    %
    % Syntax:
    %
    % Input(s):
    %
    % Output(s):
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    % See also .

    % Copyright 2023 Richard J. Cui. Created: Sun 09/17/2023 10:48:44.776 PM
    % $Revision: 0.1 $  $Date: Sun 09/17/2023 10:48:44.797 PM $
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
        in_time (1, :) double
    end % positional

    arguments
        options.InUnit (1, 1) char ...
            {mustBeMember(options.InUnit, {'index', 'uutc', 'msec', 'second', 'minute', 'hour', 'day'})} = 'uutc'
        options.OutUnit (1, 1) char ...
            {mustBeMember(options.OutUnit, {'index', 'uutc', 'msec', 'second', 'minute', 'hour', 'day'})} = 'second'
    end % optional

    in_unit = options.InUnit;
    out_unit = options.OutUnit;

    % =========================================================================
    % main
    % =========================================================================
    in_time_abs = this.relative2absTimePoint(in_time, in_unit);
    out_time_abs = this.SampleUnitConvert(in_time_abs, in_unit, out_unit); % TODO
    out_time = this.abs2relativeTimePoint(out_time_abs, out_unit);

end % function SessionUnitConvert

% [EOF]
