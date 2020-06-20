function mayo_out = read_mayo_mef30(varargin)

% READ_MAYO_MEF30 read header, event and data from the files formatted in MEF 3.0
%
% Syntax:
%   hdr = mayo_mef30(filename)
%   hdr = mayo_mef30(filename, password)
%   hdr = mayo_mef30(filename, password, sortchannel)
%   evt = mayo_mef30(filename, password, sortchannel, hdr)
%   dat = mayo_mef30(filename, password, sortchannel, hdr, begsample, endsample, chanindx)
%
% Input(s):
%   filename        - [char] name of the file or folder of the dataset
%   password        - [struct] (opt) password structure of MEF 3.0 data (see MEFSession_3p0)
%   sortchannel     - [char] (opt) sort channel order either alphabetically 'alphabet' or 
%                     numerically 'number' (default = 'alphabet')
%   hdr             - [struct] (opt) header structure of the dataset (see FT_READ_HEADER; default = struct([]))
%   begsample       - [num] (opt) first sample to read (default = [])
%   endsample       - [num] (opt) last smaple to read (default = [])
%   chanindx        - [num] (opt) list of channel indices to read (default = [])
%
% Output(s):
%   hdr             - [struct] header structure of the dataset (see FT_READ_HEADER)
%   evt             - [struct] event structure of the dataset (see FT_READ_EVENT)
%   dat             - [num] data read in
%
% Example:
%
% Note:
%
% References:
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_EVENT, FT_READ_DATA

% Copyright 2020 Richard J. Cui. Created: Sat 03/21/2020  5:26:02.846 PM
% $Revision: 0.5 $  $Date: Thu 04/09/2020 10:44:31.535 PM $
%
% Multimodel Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(varargin{:});
filename    = q.filename;
password    = q.password;
if isempty(password)
    password = struct('Level1Password', '', 'Level2Password', '',...
        'AccessLevel', 1);
end % if
sortchannel = q.sortchannel;
hdr         = q.hdr;
begsample   = q.begsample;
endsample   = q.endsample;
chanindx    = q.chanindx;

% =========================================================================
% main
% =========================================================================
% check the consistency of SortChannel
% ------------------------------------
if isempty(sortchannel)
    if isempty(hdr)
        sortchannel = 'alphabet';
    else
        sortchannel = hdr.SortChannel;
    end % if
end % if

if ~isempty(hdr) && ~strcmpi(hdr.SortChannel, sortchannel)
    warning('off', 'backtrace')
    warning('mayo_mef30:invalidSortChannel',...
        'SortChannel provided -%s- is not consistent with -%s- in header. use that in header',...
        sortchannel, hdr.SortChannel)
    warning('on', 'backtrace')
    
    sortchannel = hdr.SortChannel;
end % if

% setup the instance of the object
% --------------------------------
mef_ft = MEFFieldTrip_3p0(filename, password, 'SortChannel', sortchannel); % dealing MEF 3.0 data for FieldTrip
channames = mef_ft.SelectedChannel;

% get the desired information
% ---------------------------
switch nargin
    case {1, 2, 3}
        % get header
        mayo_out = mef_ft.getHeader(channames);
        mayo_out.SortChannel = sortchannel;
    case 4
        % get event
        mayo_out = mef_ft.getEvent(channames);
    case 7
        % get data
        mayo_out = mef_ft.getData(channames, hdr.sampleunit, begsample, endsample,...
            chanindx);
    otherwise
        % error
        error('FieldTrip:mayo_mef30:invalidInput',...
            'invalid number of inputs of the function')
end % switch

end % function mayo_mef30

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% default
default_pw = struct([]);
default_sc = 'alphabet'; % sort channel
default_hr = struct([]);
default_bs = [];
default_es = [];
default_ci = [];

% parse rule
p = inputParser;
p.addRequired('filename', @ischar);
p.addOptional('password', default_pw, @(x) isstruct(x) || isempty(x));
p.addOptional('sortchannel', default_sc, @ischar);
p.addOptional('hdr', default_hr, @isstruct);
p.addOptional('begsample', default_bs, @isnumeric);
p.addOptional('endsample', default_es, @isnumeric);
p.addOptional('chanindx', default_ci, @isnumeric);

% parse and return
p.parse(varargin{:});
q = p.Results;

end % funciton

% [EOF]