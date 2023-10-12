function out_time = SampleUnitConvert(this, in_time, varargin)
% MULTISCALEELECTROPHYSIOLOGYFILE_3P0.SAMPLEUNITCONVERT convert units of time points
% 
% Syntax:
%   out_time = SampleUnitConvert(this, in_time)
%   out_time = SampleUnitConvert(__, in_unit)
%   out_time = SampleUnitConvert(__, in_unit, out_unit)
% 
% Input(s):
%   this            - [obj] MultiscaleElectrophysiologyFile_2p1 object
%   in_time         - [arr] input time points in in_unit
%   in_unit         - [str] (opt) unit of in_time (default: uUTC)
%   out_unit        - [str] (opt) unit of out_time (default: second)
% 
% Output(s):
%   out_time        - [arr] output time points in out_unit
% 
% Note:
%   in_time and out_time are absolute time points.
% 
% See also SessionUnitConvert.

% Copyright 2020 Richard J. Cui. Created: Mon 01/20/2020 10:06:49.507 PM
% $Revision: 0.2 $  $Date: Wed 01/22/2020 10:02:38.457 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, in_time, varargin{:});
in_time = q.in_time;
in_unit = q.in_unit;
out_unit = q.out_unit;

% =========================================================================
% main
% =========================================================================
if strcmpi(in_unit, out_unit) == true
    out_time = in_time;
elseif strcmpi(in_unit, 'index') == true
    out_time = this.SampleIndex2Time(in_time, out_unit);
elseif strcmpi(out_unit, 'index') == true
    out_time = this.SampleTime2Index(in_time, in_unit);
else
    out_index = this.SampleTime2Index(in_time, in_unit);
    out_time = this.SampleIndex2Time(out_index, out_unit);
end % if

end

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% default
default_inut = 'uutc';
default_outut = 'second';
expected_ut = {'index', 'uutc', 'msec', 'second', 'minute', 'hour', 'day'};

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addRequired('in_time', @isnumeric);
p.addOptional('in_unit', default_inut,...
    @(x) any(validatestring(x, expected_ut)));
p.addOptional('out_unit', default_outut,...
    @(x) any(validatestring(x, expected_ut)));

% parse and return the results
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]