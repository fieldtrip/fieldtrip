function [x, t] = importSignal(this, varargin)
% MULTISCALEELECTROPHYSIOLOGYFILE_2P1.IMPORTMEF Import MEF channel into MATLAB
% 
% Syntax:
%   [x, t] = importSignal(this)
%   [x, t] = importSignal(__, start_end)
%   [x, t] = importSignal(__, start_end, st_unit)
%   [x, t] = importSignal(__, start_end, st_unit, filepath) 
%   [x, t] = importSignal(__, start_end, st_unit, filepath, filename)
%   [x, t] = importSignal(__, 'SubjectPassword', subj_pw)
%   [x, t] = importSignal(__, 'SessionPassword', sess_pw)
%   [x, t] = importSignal(__, 'DataPassword', data_pw)
% 
% Imput(s):
%   this            - [obj] MultiscaleElectrophysiologyFile object
%   start_end       - [1 x 2 array] (opt) [start time/index, end time/index] of 
%                     the signal to be extracted fromt the file (default:
%                     the entire signal)
%   st_unit         - [str] (opt) unit of start_end: 'Index' (default), 'uUTC',
%                     'Second', 'Minute', 'Hour', and 'Day'
%   filepath        - [str] (opt) directory of the session
%   filename        = [str] (opt) filename of the channel
%   subj_pw         - [str] (para) subject password
%   sess_pw         - [str] (para) session password
%   data_pw         - [str] (para) data password
% 
% Output(s):
%   x               - [num array] extracted signal
%   t               - [num array] time indices of the signal in the file
% 
% Note:
%   Import data from one channel.
% 
% See also .

% Copyright 2019-2020 Richard J. Cui. Created: Mon 04/29/2019 10:33:58.517 PM
% $Revision: 0.9 $  $Date: Thu 01/09/2020  4:11:13.040 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, varargin{:});
start_end = q.start_end;
st_unit = q.st_unit;
switch lower(st_unit)
    case 'index'
        se_index = start_end;
    otherwise
        se_index = this.SampleTime2Index(start_end, st_unit);
end % switch
filepath = q.filepath;
if isempty(filepath)
    filepath = this.FilePath;
else
    this.FilePath = filepath;
end % if
filename = q.filename;
if isempty(filename)
    filename = this.FileName;
else
    this.FileName = filename;
end % if

subj_pw = q.SubjectPassword;
if isempty(subj_pw) == false
    this.SubjectPassword = subj_pw;
end % if

sess_pw = q.SessionPassword;
if isempty(sess_pw) == false
    this.SessionPassword = sess_pw;
end % if

data_pw = q.DataPassword;
if isempty(data_pw) == false
    this.DataPassword = data_pw;
end % if

% check
if se_index(1) < 1
    se_index(1) = 1; 
    warning('MultiscaleElectrophysiologyFile_2p1:ImportSignal',...
        'Reqested data samples before the recording are discarded')
end % if
if se_index(2) > this.Header.number_of_samples
    se_index(2) = this.Header.number_of_samples; 
    warning('MultiscaleElectrophysiologyFile_2p1:ImportSignal',...
        'Reqested data samples after the recording are discarded')
end % if
wholename = fullfile(filepath, filename);

% verbose
num_samples = diff(start_end)+1;
if num_samples > 2^20
    verbo = true; 
else
    verbo = false;
end % if

% =========================================================================
% load the data
% =========================================================================
pw = this.SessionPassword;
if verbo, fprintf('-->Loading...'), end % if
x = decompress_mef_2p1(wholename, se_index(1), se_index(2), pw);
x = double(x(:)).'; % change to row vector
% find the indices corresponding to physically collected data
if nargout == 2
    t = se_index(1):se_index(2);
end % if
if verbo, fprintf('Done!\n'), end % if

end

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(this, varargin)

% defaults
start_ind = this.SampleTime2Index(this.Header.recording_start_time);
end_ind = this.SampleTime2Index(this.Header.recording_end_time);
defaultSE = [start_ind, end_ind];
defaultSTUnit = 'index';
expectedSTUnit = {'index', 'uutc', 'second', 'minute', 'hour', 'day'};
default_fp = '';
default_fn = '';
default_subjpw = '';
default_sesspw = '';
default_datapw = '';

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addOptional('start_end', defaultSE,...
    @(x) isnumeric(x) & numel(x) == 2 & x(1) <= x(2));
p.addOptional('st_unit', defaultSTUnit,...
    @(x) any(validatestring(x, expectedSTUnit)));
p.addOptional('filepath', default_fp, @isstr);
p.addOptional('filename', default_fn, @isstr);
p.addParameter('SubjectPassword', default_subjpw, @isstr);
p.addParameter('SessionPassword', default_sesspw, @isstr);
p.addParameter('DataPassword', default_datapw, @isstr);

% parse and return the results
p.parse(this, varargin{:});
q = p.Results;

end % function

% [EOF]