function out_time = SessionUnitConvert(this, in_time, varargin)
% MEFSESSION.SESSIONUNITCONVERT convert units of relative time points
% 
% Syntax:
%   out_time = SessionUnitConvert(this, in_time)
%   out_time = SessionUnitConvert(__, in_unit)
%   out_time = SessionUnitConvert(__, in_unit, out_unit)
% 
% Input(s):
%   this            - [obj] MEFSession_2p1 object
%   in_time         - [arr] input time points in in_unit
%   in_unit         - [str] (opt) unit of in_time (default: uUTC)
%   out_unit        - [str] (opt) unit of out_time (default: second)
% 
% Output(s):
%   out_time        - [arr] output time points in out_unit
% 
% Note:
%   in_time and out_time are relative time points.
% 
% See also SampleUnitConvert.

% Copyright 2020 Richard J. Cui. Created: Wed 01/22/2020 10:02:38.457 PM
% $Revision: 0.2 $  $Date: Sun 03/22/2020  9:15:46.097 PM $
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
in_time_abs = this.relative2absTimePoint(in_time, in_unit);
out_time_abs = this.SampleUnitConvert(in_time_abs, in_unit, out_unit);
out_time = this.abs2relativeTimePoint(out_time_abs, out_unit);

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