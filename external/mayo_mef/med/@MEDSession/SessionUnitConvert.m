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
    % $Revision: 0.2 $  $Date: Mon 10/09/2023 12:57:35.798 AM $
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
        options.InUnit (1, :) char ...
            {mustBeMember(options.InUnit, {'index', 'uutc', 'msec', 'second', 'minute', 'hour', 'day'})} = 'uutc'
        options.OutUnit (1, :) char ...
            {mustBeMember(options.OutUnit, {'index', 'uutc', 'msec', 'second', 'minute', 'hour', 'day'})} = 'second'
    end % optional

    in_unit = options.InUnit;
    out_unit = options.OutUnit;

    % =========================================================================
    % main
    % =========================================================================
    out_time = this.SampleUnitConvert(in_time, in_unit, out_unit);

end % function SessionUnitConvert

% [EOF]
