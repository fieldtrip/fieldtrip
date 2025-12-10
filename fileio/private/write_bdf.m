function write_bdf(filename, data, fs, channel_names, varargin)

% WRITE_BDF Write EEG data to a Biosemi BDF file
%
%   write_bdf(filename, data, fs, channel_names)
%   write_bdf(filename, data, fs, channel_names, 'Param', value, ...)
%
% Inputs:
%   filename      - Output BDF filename (including .bdf extension)
%   data          - EEG data matrix (Nchannels × Ntimepoints)
%   fs            - Sampling frequency (Hz)
%   channel_names - Cell array of channel names (length Nchannels)
%
% Optional parameters:
%   'SubjectID'   - Subject identification string (default: 'X')
%   'RecordingID' - Recording identification string (default: 'X')
%   'PhysicalMax' - Physical maximum for each channel (default: 32767)
%   'PhysicalMin' - Physical minimum for each channel (default: -32768)
%   'DigitalMax'  - Digital maximum for each channel (default: 8388607)
%   'DigitalMin'  - Digital minimum for each channel (default: -8388608)
%   'ScaleFactor' - Scaling factor for each channel (default: 1)
%   'Transducer'  - Transducer type for each channel (default: 'Active electrode')
%   'Prefilter'   - Prefiltering for each channel (default: 'HP:0.16Hz LP:500Hz')
%
% Example:
%   data = randn(32, 1000); % 32 channels, 1000 timepoints
%   fs = 256; % 256 Hz sampling rate
%   ch_names = arrayfun(@(x) sprintf('EEG%02d', x), 1:32, 'UniformOutput', false);
%   write_bdf('test.bdf', data, fs, ch_names, 'SubjectID', 'SUBJ01');
%
% See also READ_BDF

% This code was generated using DeepSeek and subsequently tested and fixed by Robert
% Oostenveld, hence I don't know the appropriate copyrights.
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org for the
% documentation and details.
%
% $Id$

% Parse optional parameters
p = inputParser;
addParameter(p, 'SubjectID', '?', @ischar);
addParameter(p, 'RecordingID', '?', @ischar);
addParameter(p, 'PhysicalMax', 32767, @isnumeric);
addParameter(p, 'PhysicalMin', -32768, @isnumeric);
addParameter(p, 'DigitalMax', 8388607, @isnumeric);
addParameter(p, 'DigitalMin', -8388608, @isnumeric);
addParameter(p, 'ScaleFactor', 1, @isnumeric);
addParameter(p, 'Transducer', 'Active electrode', @(x) ischar(x) || iscell(x));
addParameter(p, 'Prefilter', 'HP:?Hz LP:?Hz', @(x) ischar(x) || iscell(x));
parse(p, varargin{:});

% Get parameters
subject_id = p.Results.SubjectID;
recording_id = p.Results.RecordingID;
physical_max = p.Results.PhysicalMax;
physical_min = p.Results.PhysicalMin;
digital_max = p.Results.DigitalMax;
digital_min = p.Results.DigitalMin;
scale_factor = p.Results.ScaleFactor;
transducer = p.Results.Transducer;
prefilter = p.Results.Prefilter;

% Ensure inputs are correct
[Nchannels, Ntimepoints] = size(data);
if Nchannels ~= length(channel_names)
  error('Number of channels in data (%d) does not match length of channel_names (%d)', ...
    Nchannels, length(channel_names));
end

% Handle scalar vs. vector parameters
if isscalar(physical_max), physical_max = repmat(physical_max, Nchannels, 1); end
if isscalar(physical_min), physical_min = repmat(physical_min, Nchannels, 1); end
if isscalar(digital_max), digital_max = repmat(digital_max, Nchannels, 1); end
if isscalar(digital_min), digital_min = repmat(digital_min, Nchannels, 1); end
if isscalar(scale_factor), scale_factor = repmat(scale_factor, Nchannels, 1); end
if ischar(transducer), transducer = repmat({transducer}, Nchannels, 1); end
if ischar(prefilter), prefilter = repmat({prefilter}, Nchannels, 1); end

% Open file for writing (binary, little-endian)
fid = fopen(filename, 'w', 'ieee-le');
if fid == -1
  error('Could not open file %s for writing', filename);
end

% Calculate number of records needed
record_duration = 1; % 1-second records are common
samples_per_record = fs * record_duration;
num_records = ceil(Ntimepoints / samples_per_record);

% ====================
% Write BDF header
% ====================

% 1. Version (8 bytes)
fwrite(fid, [255 'BIOSEMI'], 'uint8'); % prevent interpretation as UTF8

% 2. Patient ID (80 bytes)
fwrite(fid, sprintf('%-80s', subject_id), 'char');

% 3. Recording ID (80 bytes)
fwrite(fid, sprintf('%-80s', recording_id), 'char');

% 4. Start date (8 bytes) - format dd.mm.yy
start_date = datestr(now, 'dd.mm.yy');
fwrite(fid, sprintf('%-8s', start_date), 'char');

% 5. Start time (8 bytes) - format hh.mm.ss
start_time = datestr(now, 'HH.MM.SS');
fwrite(fid, sprintf('%-8s', start_time), 'char');

% 6. Header size (8 bytes) - 256 + (256 * Nchannels)
header_size = 256 + (256 * Nchannels);
fwrite(fid, sprintf('%-8d', header_size), 'char');

% 7. Reserved (44 bytes)
fwrite(fid, sprintf('%-44s', ' '), 'char');

% 8. Number of data records (8 bytes) - -1 if unknown
fwrite(fid, sprintf('%-8d', num_records), 'char');

% 9. Duration of data record (8 bytes) - in seconds
fwrite(fid, sprintf('%-8d', record_duration), 'char');

% 10. Number of channels (4 bytes)
fwrite(fid, sprintf('%-4d', Nchannels), 'char');


% ====================
% Channel headers (256 bytes per channel)
% ====================

for ch = 1:Nchannels
  % 1. Channel label (16 bytes)
  fwrite(fid, sprintf('%-16s', channel_names{ch}), 'char');
end

for ch = 1:Nchannels
  % 2. Transducer type (80 bytes)
  fwrite(fid, sprintf('%-80s', transducer{ch}), 'char');
end

for ch = 1:Nchannels
  % 3. Physical dimension (uV) (8 bytes)
  fwrite(fid, sprintf('%-8s', 'uV'), 'char');
end

for ch = 1:Nchannels
  % 4. Physical minimum (8 bytes)
  fwrite(fid, sprintf('%-8g', physical_min(ch)), 'char');
end

for ch = 1:Nchannels
  % 5. Physical maximum (8 bytes)
  fwrite(fid, sprintf('%-8g', physical_max(ch)), 'char');
end

for ch = 1:Nchannels
  % 6. Digital minimum (8 bytes)
  fwrite(fid, sprintf('%-8d', digital_min(ch)), 'char');
end

for ch = 1:Nchannels
  % 7. Digital maximum (8 bytes)
  fwrite(fid, sprintf('%-8d', digital_max(ch)), 'char');
end

for ch = 1:Nchannels
  % 8. Prefiltering (80 bytes)
  fwrite(fid, sprintf('%-80s', prefilter{ch}), 'char');
end

for ch = 1:Nchannels
  % 9. Number of samples per record (8 bytes)
  fwrite(fid, sprintf('%-8d', samples_per_record), 'char');
end

for ch = 1:Nchannels
  % 10. Reserved (32 bytes)
  fwrite(fid, sprintf('%-32s', ' '), 'char');
end

% ====================
% Write data records
% ====================

% Pad data with zeros if needed
if mod(Ntimepoints, samples_per_record) ~= 0
  pad_size = num_records * samples_per_record - Ntimepoints;
  data = [data, zeros(Nchannels, pad_size)];
end

% Reshape data into records
data_records = reshape(data, Nchannels, samples_per_record, num_records);

% Write each record
for rec = 1:num_records
  for ch = 1:Nchannels
    % Scale the data to digital units
    ch_data = data_records(ch, :, rec) * scale_factor(ch);

    % Convert to 24-bit integers (3 bytes per sample)
    ch_data = int32(ch_data);

    % Write to file (3 bytes per sample, little-endian)
    n = fwrite(fid, ch_data, 'bit24');
  end
end

% Close file
fclose(fid);

fprintf('Successfully wrote %d channels of EEG data to %s\n', Nchannels, filename);
end
