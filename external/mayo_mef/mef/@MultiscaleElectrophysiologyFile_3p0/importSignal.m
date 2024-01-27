function [x, t] = importSignal(this, varargin)
% MULTISCALEELECTROPHYSIOLOGYFILE_3P0.IMPORTMEF Import MEF 3.0 channel into MATLAB
% 
% Syntax:
%   [x, t] = importSignal(this)
%   [x, t] = importSignal(__, start_end)
%   [x, t] = importSignal(__, start_end, st_unit)
%   [x, t] = importSignal(__, start_end, st_unit, filepath) 
%   [x, t] = importSignal(__, start_end, st_unit, filepath, filename)
%   [x, t] = importSignal(__, 'Level1Password', level_1_pw)
%   [x, t] = importSignal(__, 'Level2Password', level_2_pw)
%   [x, t] = importSignal(__, 'AccessLevel', access_level)
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
%   level_1_pw      - [str] (para) password of level 1 (default = this.Level1Password)
%   level_2_pw      - [str] (para) password of level 2 (default = this.Level2Password)
%   access_level    - [str] (para) data decode level to be used
%                     (default = this.AccessLevel)
% 
% Output(s):
%   x               - [num array] extracted signal
%   t               - [num array] time indices of the signal in the file
% 
% Note:
%   Import data from one channel of MEF 3.0 file into MatLab.
% 
% See also .

% Copyright 2020 Richard J. Cui. Created: Wed 02/05/2020 10:24:56.722 PM
% $Revision: 0.3 $  $Date: Fri 04/03/2020  5:11:57.375 PM $
%
% Multimodel Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
% parse inputs now
% ----------------
q = parseInputs(this, varargin{:});
start_end = q.start_end;
st_unit = q.st_unit;
filepath = q.filepath;
filename = q.filename;
l1_pw = q.Level1Password;
l2_pw = q.Level2Password;
al = q.AccessLevel;

% password
% --------
if isnan(l1_pw)
    l1_pw = this.Level1Password;
else
    this.Level1Password = l1_pw;
end % if

if isnan(l2_pw)
    l2_pw = this.Level2Password;
else
    this.Level2Password = l2_pw;
end % if

if isnan(al)
    al = this.AccessLevel;
else
    this.AccessLevel = al;
end % if

pw = this.processPassword('Level1Password', l1_pw,...
                          'Level2Password', l2_pw,...
                          'AccessLevel', al);

% get the channel metadata if both filepath and filename are provided
% -------------------------------------------------------------------
if isempty(filepath)
    filepath = this.FilePath;
else
    this.FilePath = filepath;
end % if
if isempty(filename)
    filename = this.FileName;
else
    this.FileName = filename;
end % if

wholename = fullfile(filepath, filename);
if ~isempty(filepath) && ~isempty(filename)
    [header, channel] = this.readHeader(wholename, pw, al);
    this.Header = header;
    this.Channel = channel; 
end % if

% start and end time points
% -------------------------
switch lower(st_unit)
    case 'index'
        se_index = start_end;
    otherwise
        se_index = this.SampleTime2Index(start_end, st_unit);
end % switch

if isempty(start_end) == true
    start_ind = this.SampleTime2Index(this.Channel.earliest_start_time);
    end_ind = this.SampleTime2Index(this.Channel.latest_end_time);
    se_index = [start_ind, end_ind];
end % if

% check
if se_index(1) < 1
    se_index(1) = 1; 
    warning('MultiscaleElectrophysiologyFile_3p0:ImportSignal:discardSample',...
        'Reqested data samples before the recording are discarded')
end % if
if se_index(2) > this.Channel.metadata.section_2.number_of_samples
    se_index(2) = this.Channel.metadata.section_2.number_of_samples; 
    warning('MultiscaleElectrophysiologyFile_3p0:ImportSignal:discardSample',...
        'Reqested data samples after the recording are discarded')
end % if

% verbose
% -------
num_samples = diff(start_end)+1;
if num_samples > 2^20
    verbo = true; 
else
    verbo = false;
end % if

% =========================================================================
% load the data
% =========================================================================
if verbo
    [~, thisChannel] = fileparts(wholename);
    fprintf(['-->Loading ' thisChannel ' ...'])
    clear thisChannel
end % if
x = this.read_mef_ts_data_3p0(wholename, pw, 'samples', se_index(1), se_index(2));
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
defaultSE = [];
defaultSTUnit = 'index';
expectedSTUnit = {'index', 'uutc', 'second', 'minute', 'hour', 'day'};
default_fp = '';
default_fn = '';
default_l1pw = NaN;
default_l2pw = NaN;
default_al = nan; % access level

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addOptional('start_end', defaultSE,...
    @(x) isnumeric(x) & numel(x) == 2 & x(1) <= x(2));
p.addOptional('st_unit', defaultSTUnit,...
    @(x) any(validatestring(x, expectedSTUnit)));
p.addOptional('filepath', default_fp, @isstr);
p.addOptional('filename', default_fn, @isstr);
p.addParameter('Level1Password', default_l1pw, @(x) ischar(x) || isnan(x));
p.addParameter('Level2Password', default_l2pw, @(x) ischar(x) || isnan(x));
p.addParameter('AccessLevel', default_al, @isnumeric);

% parse and return the results
p.parse(this, varargin{:});
q = p.Results;

end % function

% [EOF]