function rel_time = abs2relativeTimePoint(this, abs_time, unit)
% MEFSESSION.ABS2RELATIVETIMEPOINT convert absolute time point to relative one
% 
% Syntax:
%   rel_time = abs2relativeTimePoint(this, abs_time, unit)
% 
% Input(s):
%   this            - [obj] MEFSession_2p1 object
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

% Copyright 2020 Richard J. Cui. Created: Mon 01/20/2020  4:30:22.035 PM
% $Revision: 0.1 $  $Date: Mon 01/20/2020  4:30:22.035 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, abs_time, unit);
abs_time = q.abs_time;
unit = q.unit;

% =========================================================================
% main
% =========================================================================
offset = this.getSessionRecordOffset(unit);
rel_time = abs_time-offset;

end

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% default
expected_ut = {'index', 'uutc', 'msec', 'second', 'minute', 'hour', 'day'};

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addRequired('abs_time', @isnumeric);
p.addRequired('unit',  @(x) any(validatestring(x, expected_ut)));

% parse and return the results
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]