function [header, channel] = readHeader(this, varargin)
% MULTISCALEELECTROPHYSIOLOGYFILE_3P0.READHEADER Read UNIVERSAL HEADER structure from MEF 3.0
% 
% Syntax:
%   [header, channel] = readHeader(this)
%   [header, channel] = readHeader(this, wholename)
%   [header, channel] = readHeader(this, wholename, password)
%   [header, channel] = readHeader(this, wholename, password, access_level)
% 
% Input(s):
%   this            - [obj] MultiscaleElectrophysiologyFile_3p0 object
%   wholename       - [str] (opt) filepath + filename of MEF channel file
%   password        - [str] (opt) password of the data
%   access_level    - [num] (opt) access level of data (default = 1)
% 
% Output(s):
%   header          - [struct] UNIVERSAL_HEADER info of MEF 3.0 file
%                     .header_CRC
%                     .body_CRC
%                     .file_type_string
%                     .mef_version_major
%                     .mef_version_minor
%                     .byte_order_code
%                     .start_time
%                     .end_time
%                     .number_of_entries
%                     .maximum_entry_size
%                     .segment_number
%                     .channel_name
%                     .session_name
%                     .anonymized_name
%                     .level_UUID
%                     .file_UUID
%                     .provenance_UUID
%                     .level_1_password_validation_field
%                     .level_2_password_validation_field
%                     .protected_region
%                     .discretionary_region
%   channel         - [struct] CHAANEL info of MEF3.0 file
%                     .channel_type
%                     .metadata
%                     .number_of_segments
%                     .segments
%                     .path
%                     .name
%                     .extension
%                     .session_name
%                     .level_UUID
%                     .anonymized_name
%                     .maximum_number_of_records
%                     .maximum_record_bytes
%                     .earliest_start_time
%                     .latest_end_time
%                     .records
% 
% Note:
%   See the details of MEF file at https://github.com/benbrinkmann/mef_lib_2_1
% 
% See also read_mef_header_mex_3p0.m.

% Copyright 2020 Richard J. Cui. Created: Tue 02/04/2020  3:33:28.609 PM
% $Revision: 0.4 $  $Date: Fri 06/05/2020 10:37:43.416 AM $
%
% Multimodel Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905, USA
%
% Email: richard.cui@utoronto.ca

q = parseInputs(varargin{:});

if isempty(q)
    wholename = fullfile(this.FilePath, this.FileName);
    pw = this.Password;
else
    wholename = fullfile(q.filepath, q.filename);
    pw = q.password;
    al = q.access_level;
    % update
    this.FilePath = q.filepath;
    this.FileName = q.filename;
    switch al
        case 1
            this.Level1Password = pw;
        case 2
            this.Level2Password = pw;
    end % switch
    this.AccessLevel = al;
end % if

% get the channel info
map_tsi = true;
[sess_path,chan_name] = fileparts(wholename);
if isempty(pw)
    md = read_mef_session_metadata(sess_path,[],map_tsi);
else
    md = read_mef_session_metadata(sess_path,pw,map_tsi);
end % if
all_chan_name = {md.time_series_channels.name};
chan_indx = ismember(all_chan_name,chan_name);

% get the header and channel
channel = md.time_series_channels(chan_indx);
header  = channel.segments(1).time_series_data_uh;

end %function

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% defaults
default_pw = '';
default_al = 1;

% parse rules
p = inputParser;
p.addOptional('wholename', '', @isstr);
p.addOptional('password', default_pw, @isstr);
p.addOptional('access_level', default_al, @isnumeric);

% parse and return the results
p.parse(varargin{:});
if isempty(p.Results.wholename)
    q = [];
else
    [fp, fn, ext] = fileparts(p.Results.wholename);
    q.filepath = fp;
    q.filename = [fn, ext];
    q.password = p.Results.password;
    q.access_level = p.Results.access_level;
end % if

end % function

% [EOF]

