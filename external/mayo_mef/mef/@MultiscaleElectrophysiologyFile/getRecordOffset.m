function record_offset = getRecordOffset(this, unit)
% MULTISCALEELECTROPHYSIOLOGY.GETRECORDOFFSET get offset time of recording in specified unit
% 
% Syntax:
%   record_offset = getRecordOffset(this, unit)
% 
% Input(s):
%   this            - [obj] MEFSession_2p1 object
%   unit            - [char] unit of the output offset
% 
% Output(s):
%   record_offset   - [num] recording offset of time
% 
% Note:
% 
% Seaa slso .

% Copyright 2020 Richard J. Cui. Created: Mon 01/20/2020  4:30:22.035 PM
% $Revision: 0.2 $  $Date: Thu 02/06/2020  9:23:51.173 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, unit);
unit = q.unit;

% =========================================================================
% main
% =========================================================================
if strcmpi(unit, 'index') % get recoridng start time in unit
    record_offset = 0;
else
    record_offset = this.SampleIndex2Time(1, unit);
end % if

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
p.addRequired('unit',  @(x) any(validatestring(x, expected_ut)));

% parse and return the results
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]