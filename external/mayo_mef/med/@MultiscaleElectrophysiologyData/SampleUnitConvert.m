function out_time = SampleUnitConvert(this, in_time, in_unit, out_unit)
    % MULTISCALEELECTROPHYSIOLOGYDATA.SAMPLEUNITCONVERT convert units of time stamps
    %
    % Syntax:
    %   out_time = SampleUnitConvert(this, in_time)
    %   out_time = SampleUnitConvert(__, in_unit)
    %   out_time = SampleUnitConvert(__, in_unit, out_unit)
    %
    % Input(s):
    %   this            - [obj] MultiscaleElectrophysiologyData object
    %   in_time         - [arr] input time points in in_unit
    %   in_unit         - [str] (opt) unit of in_time (default: uutc)
    %   out_unit        - [str] (opt) unit of out_time (default: second)
    %
    % Output(s):
    %   out_time        - [arr] output time points in out_unit
    %
    % Note:
    %   in_time and out_time are absolute time points.
    %
    % See also SessionUnitConvert.

    % Copyright 2023 Richard J. Cui. Created: Mon 10/09/2023  9:46:29.009 PM
    % $Revision: 0.1 $  $Date: Mon 10/09/2023  9:46:29.014 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        this (1, 1) MultiscaleElectrophysiologyData
        in_time (1, :) double
        in_unit (1, 1) string = "uutc"
        out_unit (1, 1) string = "second"
    end % positional

    % ======================================================================
    % main
    % ======================================================================
    if strcmpi(in_unit, out_unit) == true
        out_time = in_time;
    elseif strcmpi(in_unit, 'index') == true
        out_time = this.SampleIndex2Time(in_time, st_unit = out_unit);
    elseif strcmpi(out_unit, 'index') == true
        out_time = this.SampleTime2Index(in_time, st_unit = in_unit);
    else
        out_index = this.SampleTime2Index(in_time, st_unit = in_unit);
        out_time = this.SampleIndex2Time(out_index, st_unit = out_unit);
    end % if

end % function SampleUnitConvert

% [EOF]
