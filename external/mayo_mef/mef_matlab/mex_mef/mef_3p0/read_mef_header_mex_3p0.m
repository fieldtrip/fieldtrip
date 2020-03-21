function [header, channel] = read_mef_header_mex_3p0(channel_path,...
    password, map_indices_flag)
% READ_MEF_HEADER_MEX_3P0 read universal header and channel metadata
% 
% Syntax:
%   [header, channel] = read_mef_header_mex_3p0(channel_path)
%   [header, channel] = read_mef_header_mex_3p0(__, password)
%   [header, channel] = read_mef_header_mex_3p0(__, password, map_indices_flag)
% 
% Input(s):
%   channel_path    - [str] character vector of channel path
%   password        - [str] (opt) character string of password (defulat =
%                     '')
%   map_indices_flag- [num] (opt) (default = 1)
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
% See also .

% Copyright 2020 Richard J. Cui. Created:Tue 02/04/2020 12:04:38.033 PM
% $Revision: 0.1 $  $Date: Tue 02/04/2020 12:04:38.033 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% message
% =========================================================================
warning(['You probably haven''t built the mex file for %s.c.',...
    'Please run the make file to built the binary.'], mfilename)

% =========================================================================
% dummy code
% =========================================================================
x = 1;
switch x
    case 1
        y = channel_path;
    case 2
        y = password;
    case 3
        y = map_indices_flag;
end % switch
header = y;
channel = [];

end % function

% [EOF]