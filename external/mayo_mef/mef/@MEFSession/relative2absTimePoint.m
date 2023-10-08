function abs_time = relative2absTimePoint(this, rel_time, unit)
% MEFSESSION.RELATIVE2ABSTIMEPOINT convert relative time point to absolute one
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
% Note:
% 
% Seaa slso .

% Copyright 2020 Richard J. Cui. Created: Mon 01/20/2020  4:30:22.035 PM
% $Revision: 0.2 $  $Date: Sat 02/08/2020 11:40:56.279 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, rel_time, unit);
rel_time = q.rel_time;
unit = q.unit;

% =========================================================================
% main
% =========================================================================
offset = this.getSessionRecordOffset(unit);
abs_time = rel_time+offset;

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
p.addRequired('rel_time', @isnumeric);
p.addRequired('unit',  @(x) any(validatestring(x, expected_ut)));

% parse and return the results
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]