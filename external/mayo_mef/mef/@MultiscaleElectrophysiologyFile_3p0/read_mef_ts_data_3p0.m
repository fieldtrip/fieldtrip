function data = read_mef_ts_data_3p0(this, channel_path, varargin)
% MultiscaleElectrophysiologyFile_3p0.READ_MEF_TS_DATA_3P0 Read the MEF 3.0 data from a time-series channel
%	
% Syntax:
%   data = read_mef_ts_data_3p0(this, channel_path)
%   data = read_mef_ts_data_3p0(__, password)
%   data = read_mef_ts_data_3p0(__, password, range_type)
%   data = read_mef_ts_data_3p0(__, password, range_type, begin, stop)
% 
% Input(s):
%   this            - [obj] MultiscaleElectrophysiologyFile_3p0 object
%   channel_path    - [char] (opt) path (absolute or relative) to the MEF3 
%                     channel folder
%   password        - [str] (opt) password to the MEF 3.0 data; Pass 
%                     empty string if not encrypted (default = '')
%                     .
%   range_type      - [char] (opt) modality that is used to define the 
%                     data-range to read, either 'time' or 'samples'
%                     (default = samples)
%   begin           - [array] (opt) Start-point for the reading of data 
%                     (either as a timepoint or samplenumber); Pass -1 to
%                     start at the first sample of the timeseries. if
%                     range_type is 'smaples', the first sample begins at 1
%                     using Matlab convention. (default = -1)
%   stop            - [array] (opt) End-point to stop the reading data
%                     (either as a timepoint or samplenumber); Pass -1 as
%                     value to end at the last sample of the timeseries
%                     (default = -1)
%
% Output(s): 
%   data            - A vector of doubles holding the channel data
%
% Note:
%   When the rangeType is set to 'samples', the function simply returns the
%   samples as they are found (consecutively) in the datafile, without any
%   regard for time or data gaps; Meaning that, if there is a time-gap
%   between samples, then these will not appear in the result returned. In
%   contrast, the 'time' rangeType will return the data with NaN values in
%   place for the missing samples.
%
% See also .

%   Copyright 2020, Max van den Boom (Multimodal Neuroimaging Lab, Mayo
%   Clinic, Rochester MN) <https://github.com/MaxvandenBoom/matmef>.
%   Adapted from PyMef (by Jan Cimbalnik, Matt Stead, Ben Brinkmann, and
%   Dan Crepeau)
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. This program is distributed in the hope that
%   it will be useful, but WITHOUT ANY WARRANTY; without even the implied
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
%   the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <https://www.gnu.org/licenses/>.

% Copyright 2020 Richard J. Cui. Adapted: Fri 01/31/2020 11:59:20.073 PM
% $Revision: 0.3 $  $Date: Wed 02/05/2020 11:36:04.992 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, channel_path, varargin{:});
ch_path = q.channel_path;
pw = q.password;
rtype = q.range_type;
begin = q.begin;
stop = q.stop;

if isempty(ch_path)
    ch_path = fullfile(this.FilePath, this.FileName);
end % if

if isempty(pw)
    pw = this.processPassword;
end % if

if strcmpi(rtype, 'samples') == true && begin ~= -1
    begin = begin-1; % change to python convention
end % if

% =========================================================================
% main
% =========================================================================
data = read_mef_ts_data(ch_path, pw, rtype, begin, stop);

end

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% defaults
default_cp = '';
default_pw = ''; % password
default_rt = 'samples'; % range_type
default_bg = -1; % begin
default_sp = -1; % stop

expected_type = {'samples', 'time'};

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addOptional('channel_path', default_cp, @isstr);
p.addOptional('password', default_pw, @isstr);
p.addOptional('range_type', default_rt, @(x) any(validatestring(x, expected_type)));
p.addOptional('begin', default_bg, @isnumeric);
p.addOptional('stop', default_sp, @isnumeric);

% parse and return the results
p.parse(varargin{:});
q = p.Results;

end % funciton

% [EOF]