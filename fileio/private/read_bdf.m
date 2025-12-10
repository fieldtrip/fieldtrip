function [data, header] = read_bdf(filename, varargin)
% READ_BDF Read EEG data from a Biosemi BDF file
%
%   [data, header] = read_bdf(filename)
%   [data, header] = read_bdf(filename, 'Param', value, ...)
%
% Inputs:
%   filename - BDF filename (including .bdf extension)
%
% Outputs:
%   data    - EEG data matrix (Nchannels × Ntimepoints)
%   header  - Structure containing header information
%
% Optional parameters:
%   'Channels'    - Cell array of channel names to read (default: all)
%   'TimeRange'   - [start end] in seconds (default: all)
%   'Verbose'     - Display progress info (default: true)
%
% Example:
%   [eeg, hdr] = read_bdf('eeg_data.bdf');
%   [eeg, hdr] = read_bdf('eeg_data.bdf', 'Channels', {'Fz', 'Cz', 'Pz'}, 'TimeRange', [10 20]);
%
% See also WRITE_BDF

% This code was generated using DeepSeek and subsequently tested and fixed by Robert
% Oostenveld, hence I don't know the appropriate copyrights.
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org for the
% documentation and details.
%
% $Id$

% Parse optional parameters
p = inputParser;
addParameter(p, 'Channels', 'all', @(x) iscell(x) || strcmpi(x, 'all'));
addParameter(p, 'TimeRange', 'all', @(x) isequal(x, 'all') || (isnumeric(x) && numel(x)==2));
addParameter(p, 'Verbose', true, @islogical);
parse(p, varargin{:});

% Open file for reading (binary, little-endian)
fid = fopen(filename, 'r', 'ieee-le');
if fid == -1
  error('Could not open file %s for reading', filename);
end

% Initialize header structure
header = struct();

% ====================
% Read BDF header
% ====================

% 1. Version (8 bytes)
header.version = deblank(char(fread(fid, 8, 'char')'));

% 2. Patient ID (80 bytes)
header.patient_id = deblank(char(fread(fid, 80, 'char')'));

% 3. Recording ID (80 bytes)
header.recording_id = deblank(char(fread(fid, 80, 'char')'));

% 4. Start date (8 bytes) - format dd.mm.yy
header.start_date = deblank(char(fread(fid, 8, 'char')'));

% 5. Start time (8 bytes) - format hh.mm.ss
header.start_time = deblank(char(fread(fid, 8, 'char')'));

% 6. Header size (8 bytes)
header.header_size = str2double(deblank(char(fread(fid, 8, 'char')')));

% 7. Reserved (44 bytes)
header.reserved = deblank(char(fread(fid, 44, 'char')'));

% 8. Number of data records (8 bytes)
header.num_data_records = str2double(deblank(char(fread(fid, 8, 'char')')));

% 9. Duration of data record (8 bytes) - in seconds
header.record_duration = str2double(deblank(char(fread(fid, 8, 'char')')));

% 10. Number of channels (4 bytes)
header.num_channels = str2double(deblank(char(fread(fid, 4, 'char')')));

% Initialize channel info
header.channels = struct(...
  'labels', cell(header.num_channels, 1), ...
  'transducer', cell(header.num_channels, 1), ...
  'physical_dim', cell(header.num_channels, 1), ...
  'physical_min', zeros(header.num_channels, 1), ...
  'physical_max', zeros(header.num_channels, 1), ...
  'digital_min', zeros(header.num_channels, 1), ...
  'digital_max', zeros(header.num_channels, 1), ...
  'prefiltering', cell(header.num_channels, 1), ...
  'samples_per_record', zeros(header.num_channels, 1));

% ====================
% Read channel headers
% ====================

for ch = 1:header.num_channels
  % 1. Channel label (16 bytes)
  header.channels(ch).labels = deblank(char(fread(fid, 16, 'char')'));
end

for ch = 1:header.num_channels
  % 2. Transducer type (80 bytes)
  header.channels(ch).transducer = deblank(char(fread(fid, 80, 'char')'));
end

for ch = 1:header.num_channels
  % 3. Physical dimension (8 bytes)
  header.channels(ch).physical_dim = deblank(char(fread(fid, 8, 'char')'));
end

for ch = 1:header.num_channels
  % 4. Physical minimum (8 bytes)
  header.channels(ch).physical_min = str2double(deblank(char(fread(fid, 8, 'char')')));
end

for ch = 1:header.num_channels
  % 5. Physical maximum (8 bytes)
  header.channels(ch).physical_max = str2double(deblank(char(fread(fid, 8, 'char')')));
end

for ch = 1:header.num_channels
  % 6. Digital minimum (8 bytes)
  header.channels(ch).digital_min = str2double(deblank(char(fread(fid, 8, 'char')')));
end

for ch = 1:header.num_channels
  % 7. Digital maximum (8 bytes)
  header.channels(ch).digital_max = str2double(deblank(char(fread(fid, 8, 'char')')));
end

for ch = 1:header.num_channels
  % 8. Prefiltering (80 bytes)
  header.channels(ch).prefiltering = deblank(char(fread(fid, 80, 'char')'));
end

for ch = 1:header.num_channels
  % 9. Number of samples per record (8 bytes)
  header.channels(ch).samples_per_record = str2double(deblank(char(fread(fid, 8, 'char')')));
end

for ch = 1:header.num_channels
  % 10. Reserved (32 bytes)
  fread(fid, 32, 'char');
end

if any([header.channels.samples_per_record]~=header.channels(1).samples_per_record)
  error('channels with different sampling rates are not supported')
end

% Calculate sampling rate (assuming same for all channels)
header.sample_rate = header.channels(1).samples_per_record / header.record_duration;

% Determine which channels to read
if isequal(p.Results.Channels, 'all')
  channels_to_read = 1:header.num_channels;
else
  [dum, channels_to_read] = ismember(p.Results.Channels, {header.channels.labels});
  channels_to_read = channels_to_read(channels_to_read > 0);
  if ~isempty(p.Results.Channels) && isempty(channels_to_read)
    error('None of the specified channels were found in the file');
  end
end

% Determine time range to read
if isequal(p.Results.TimeRange, 'all')
  start_record = 1;
  end_record = header.num_data_records;
else
  start_record = floor(p.Results.TimeRange(1)/header.record_duration + 1);
  end_record = ceil(p.Results.TimeRange(2)/header.record_duration);
  start_record = max(1, start_record);
  end_record = min(header.num_data_records, end_record);
end
num_records_to_read = end_record - start_record + 1;

% Calculate samples per channel per record
samples_per_record = header.channels(1).samples_per_record;
total_samples = samples_per_record * num_records_to_read;

% Initialize data matrix
data = zeros(length(channels_to_read), total_samples);

% ====================
% Read data records
% ====================

if p.Results.Verbose
  fprintf('Reading %d records (%.1f - %.1f sec) from %d channels...\n', ...
    num_records_to_read, ...
    (start_record-1)*header.record_duration, ...
    end_record*header.record_duration, ...
    length(channels_to_read));
end

% Skip to start record
bytes_per_record = sum([header.channels.samples_per_record]) * 3; % 3 bytes per sample
fseek(fid, header.header_size + (start_record-1)*bytes_per_record, 'bof');

% Read each record
for rec = 1:num_records_to_read
  if p.Results.Verbose && mod(rec, 50) == 0
    fprintf('Reading record %d of %d...\n', rec, num_records_to_read);
  end

  % Read all channels for this record (faster than channel by channel)
  all_data = fread(fid, [3, sum([header.channels.samples_per_record])], 'uint8')';

  % Convert to 24-bit integers
  for ch_idx = 1:length(channels_to_read)
    ch = channels_to_read(ch_idx);
    n_samples = header.channels(ch).samples_per_record;

    % Get the bytes for this channel
    offset = (ch-1)*n_samples;
    ch_bytes = all_data(offset + (1:n_samples), :);

    % Convert to 24-bit signed integer
    ch_data = ch_bytes(:,1) + 256*ch_bytes(:,2) + 65536*ch_bytes(:,3);
    ch_data(ch_data >= 2^23) = ch_data(ch_data >= 2^23) - 2^24; % handle signed

    % Scale to physical units
    phys_range = header.channels(ch).physical_max - header.channels(ch).physical_min;
    dig_range = header.channels(ch).digital_max - header.channels(ch).digital_min;
    scale_factor = phys_range / dig_range;

    ch_data = header.channels(ch).physical_min + (ch_data - header.channels(ch).digital_min) * scale_factor;

    % Store in output matrix
    start_sample = (rec-1)*samples_per_record + 1;
    end_sample = rec*samples_per_record;
    data(ch_idx, start_sample:end_sample) = ch_data';
  end
end

% Close file
fclose(fid);

% Trim to requested time range if needed
if ~isequal(p.Results.TimeRange, 'all')
  start_sample = floor((p.Results.TimeRange(1) - (start_record-1)*header.record_duration) * header.sample_rate) + 1;
  end_sample = ceil((p.Results.TimeRange(2) - (start_record-1)*header.record_duration) * header.sample_rate);
  start_sample = max(1, start_sample);
  end_sample = min(size(data,2), end_sample);
  data = data(:, start_sample:end_sample);
end

if p.Results.Verbose
  fprintf('Done. Read %d channels with %d samples each (%.1f sec at %.1f Hz).\n', ...
    size(data,1), size(data,2), size(data,2)/header.sample_rate, header.sample_rate);
end
end
