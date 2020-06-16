function sti = getSampleTimeInterval(this, varargin)
% MULTISCALEELECTROPHYSIOLOGYFILE.GETSAMPLETIMEINTERVAL get the boundary of sampling interval of samples
% 
% Syntax:
%   sti = getSampleTimeInterval(this)
%   sti = getSampleTimeInterval(__, fs)
% 
% Input(s):
%   this            - [obj] MultiscaleElectrophysiologyFile object
%   fs              - [num] (opt) sampling frequency (Hz)
% 
% Output(s):
%   sti             - [array] sample time inteveral = [lower, upper] (uUTC)
% 
% See also .

% Copyright 2020 Richard J. Cui. Created: Wed 02/05/2020 12:11:43.082 PM
% $Revision: 0.1 $  $Date: Wed 02/05/2020 12:11:43.082 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca


% =========================================================================
% parse the inputs
% =========================================================================
q = parseInputs(this, varargin{:});
fs = q.fs;
if isempty(fs), fs = this.ChanSamplingFreq; end % if

% =========================================================================
% main
% =========================================================================
mps = this.MPS; % microseconds per second
dt = mps/fs;
lower = floor(dt);
upper = lower + 1;
sti = [lower, upper];

% update
this.SampleTimeInterval = sti;

end

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% defaults
default_fs = [];

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addOptional('fs', default_fs, @isnumeric);

% parse and return the results
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]