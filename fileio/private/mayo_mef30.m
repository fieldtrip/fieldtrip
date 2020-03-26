function mayo_out = mayo_mef30(varargin)
% MAYO_MEF30 read header, event and data from the files formatted in MEF 3.0
%
% Syntax:
%   hdr = mayo_mef30(filename)
%   hdr = mayo_mef30(filename, password)
%   evt = mayo_mef30(filename, password, hdr)
%   dat = mayo_mef30(filename, password, hdr, begsample, endsample, chanindx)
%
% Input(s):
%   filename        - [char] name of the file or folder of the dataset
%   password        - [struct] (opt) password structure of MEF 3.0 data (see
%                     MEFSession_3p0)
%   hdr             - [struct] (opt) header structure of the dataset (see
%                     ft_read_header; default = struct([]))
%   begsample       - [num] (opt) first sample to read (default = [])
%   endsample       - [num] (opt) last smaple to read (default = [])
%   chanindx        - [num] (opt) list of channel indices to read (default
%                     = [])
%
% Output(s):
%   hdr             - [struct] header structure of the dataset (see 
%                     FT_READ_HEADER)
%   evt             - [struct] event structure of the dataset (see
%                     FT_READ_EVENT)
%   dat             - [num] data read in
%
% Example:
%
% Note:
%
% References:
%
% See also ft_filetype, ft_read_header, ft_read_event, ft_read_data.

% Copyright 2020 Richard J. Cui. Created: Sat 03/21/2020  5:26:02.846 PM
% $Revision: 0.3 $  $Date: Wed 03/25/2020  8:53:17.504 PM $
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
hdr         = q.hdr;
begsample   = q.begsample;
endsample   = q.endsample;
chanindx    = q.chanindx;

% =========================================================================
% main
% =========================================================================
% setup the instance of the object
% --------------------------------
mef_ft = MEFFieldTrip_3p0; % dealing MEF 3.0 data for FieldTrip
[sesspath, channames] = mef_ft.findSessPath(filename);
mef_ft.setSessionInfo(sesspath, password);

% get the desired information
% ---------------------------
switch nargin
    case {1, 2}
        % get header
        mayo_out = mef_ft.getHeader(channames);
    case 3
        % get event
        mayo_out = mef_ft.getEvent(channames);
    case 6
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
default_hr = struct([]);
default_bs = [];
default_es = [];
default_ci = [];

% parse rule
p = inputParser;
p.addRequired('filename', @ischar);
p.addOptional('password', default_pw, @(x) isstruct(x) || isempty(x));
p.addOptional('hdr', default_hr, @isstruct);
p.addOptional('begsample', default_bs, @isnumeric);
p.addOptional('endsample', default_es, @isnumeric);
p.addOptional('chanindx', default_ci, @isnumeric);

% parse and return
p.parse(varargin{:});
q = p.Results;

end % funciton

% [EOF]
