function header = readHeader(this, varargin)
% MULTISCALEELECTROPHYSIOLOGYFILE_2P1.READHEADER Read HEADER structure from MEF 2.1 file
% 
% Syntax:
%   header = readHeader(this)
%   header = readHeader(this, wholename)
%   header = readHeader(this, wholename, password)
% 
% Input(s):
%   this            - [obj] MultiscaleElectrophysiologyFile object
%   wholename       - [str] filepath + filename of MEF file
%   password        - [str] password of the data
% 
% Output(s):
%   header          - [struc] HEADER information of MEF 2.1 file
%                     .institution      : [$(63)] Institution
%                     .unencrypted_text_field
%                                       : [$(63)] Unencrypted Text Field (general
%                                         use)
%                     .encryption_algorithm
%                                       : [$(31)] Encryption Algorithm
%                                         "128-bit AES'
%                     .subject_encryption_used
%                                       : [ui1] Subject Encryption Use, 1
%                                         if subject encryption used, 0 if
%                                         not
%                     .session_encryption_used
%                                       : [ui1] Session Encryption Used (1
%                                         used, 0 not)
%                     .data_encryption_used
%                                       : [ui1] Data Encryption Used (1
%                                         used, 0 not)
%                     .byte_order_code  : [ui1] Byte Order Code
%                                         ('little-endian' or 'big-endian)
%                     .header_version_major
%                                       : [ui1] Header Major Version
%                     .header_version_minor
%                                       : [ui1] Header Minor Version
%                     .header_length    : [ui2] Header Length (in bytes)
%                     .session_unique_ID: [ui1] Session Unique Identifier
%                     .subject_first_name
%                                       : [$(31)] Subject First Name
%                                         ('none' if not entered)
%                     .subject_second_name
%                                       : [$(31)] Subject Middle Name
%                     .subject_third_name
%                                       : [$(31)] Subject Last Name
%                     .subject_id       : [$(31)] Subject ID
%                     .session_password : [$(15)] Session Password (15
%                                         character limit)
%                     (not read)        : [ui1] Subject Password Validation
%                                         Field
%                     (not read)        : [] Protected Region
%                                         (discretionary)
%                     (not read)        : [ui1] Session Password Validation
%                                         Field
%                     .number_of_samples: [ui8] Number of Entries (total
%                                         recorded samples in file)
%                     .channel_name     : [$(31)] Channel Name
%                     .recording_start_time
%                                       : [ui8] Recording Start Time (in
%                                         uUTC)
%                     .recording_end_time
%                                       : [ui8] Recording Eng Time (in
%                                         uUTC)
%                     .sampling_frequency
%                                       : [sf8] Sampling Frequency (-1 no
%                                         entry)
%                     .low_frequency_filter_setting
%                                       : [sf8] Low Frequency Filter
%                                         Setting (high-pass filter
%                                         setting)
%                     .high_frequency_filter_setting
%                                       : [sf8] High Frequency Filter
%                                         Setting (low-pass filter setting)
%                     .notch_filter_frequency
%                                       : [sf8] Notch Filter Frequency
%                     .voltage_conversion_factor
%                                       : [sf8] Voltage Coversion Factor
%                                         (microvolts per sample unit; 0 no
%                                         entry; ngative values inverted
%                                         voltage values)
%                     .acquisition_system
%                                       : [$(31)] Acquisition System
%                     .channel_comments : [$(127)] Channel Comments
%                     .study_comments   : [$(127)] Study Comments
%                     .physical_channel_number
%                                       : [si4] Physical Channel Number
%                                         (during acquisition)
%                     .compression_algorithm
%                                       : [$(31)] Compression Algorithm
%                                         (RED 1.0)
%                     .maximum_compressed_block_size
%                                       : [ui4] Maximum Compressed Block
%                                         Size (maximum bytes in compressed
%                                         block, including block header)
%                     .maximum_block_length
%                                       : [ui8] Maximum Block Length
%                                         (maximum number of samples in a
%                                         decompressed block)
%                     .block_interval   : [ui8] Block Interval
%                                         (microseconds between blocks; 0
%                                         variable block intervals)
%                     .maximum_data_value
%                                       : [si4] Maximum Data Value (the
%                                         largest data value in the file)
%                     .minimum_data_value
%                                       : [si4] Minimum Data Value (the
%                                         smallest data value in the file)
%                     .index_data_offset: [ui8] Offset to Block Indices
%                                         Data (offset to start of block
%                                         indices)
%                     .number_of_index_entries
%                                       : [ui8] Number of Block Index
%                                         Entries (total number data blocks
%                                         indexed; each block is indexed by
%                                         three numbers (triplets))
%                     .block_header_length
%                                       : [ui2] length of encoded data
%                                         block header in bytes
%                     .GMT_offset       : [sf4] GMT offset (GMT offset time
%                                         of file recording)
%                     .discontinuity_data_offset
%                                       : [ui8] Offset to Discontinuity
%                                         Indices Data (offset to start of
%                                         discontinuity indices which
%                                         contains block indicies where
%                                         discontinuity occred; 1st block
%                                         is numbered as zero)
%                     .number_of_discontinuity_entries
%                                       : [ui8] Number of Discontinuity
%                                         Index Entries (1 indicates the
%                                         1st block is discontinuous)
%                     (not read)        : [ui1] Unused (random bytes)
%                     .file_unique_ID   : [ui1] File Unique Identifier
%                     .anonymized_subject_name
%                                       : [$(63)] Anonymized Subject Name
%                     .header_crc       : [ui4] Header CRC (cyclically
%                                         redundant checksum)
% 
% Note:
%   See the details of MEF file at https://github.com/benbrinkmann/mef_lib_2_1
% 
% See also .

% Copyright 2019-2020 Richard J. Cui. Created: Mon 04/29/2019 10:33:58.517 PM
% $Revision: 0.7 $  $Date: Tue 02/04/2020 12:04:38.033 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

q = parseInputs(varargin{:});

if isempty(q)
    wholename = fullfile(this.FilePath, this.FileName);
    pw = this.Password;
else
    wholename = fullfile(q.filepath, q.filename);
    pw = q.password;
    % TODO: update
    this.FilePath = q.filepath;
    this.FileName = q.filename;
    this.SubjectPassword = q.password; % this is not correct
end % if

header = read_mef_header_2p1(wholename, pw);
this.Header = header;

end %function

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% defaults
default_pw = '';

% parse rules
p = inputParser;
p.addOptional('wholename', '', @isstr);
p.addOptional('password', default_pw, @isstr);

% parse and return the results
p.parse(varargin{:});
if isempty(p.Results.wholename)
    q = [];
else
    [fp, fn, ext] = fileparts(p.Results.wholename);
    q.filepath = fp;
    q.filename = [fn, ext];
    q.password = p.Results.password;
end % if

end % function

% [EOF]

