function varargout = physionet_hea(filename, hdr, begsample, endsample, chanindx)

% PHYSIONET_HEA reads data from PhysioNet .hea and .dat files
%
% Use as
%   hdr = motion_c3d(filename);
%   dat = motion_c3d(filename, hdr, begsample, endsample, chanindx);
%   evt = motion_c3d(filename, hdr);
%
% This is for data hosted on https://physionet.org following the format 
% documentd on https://wfdb.io
%
% The format specification is not very clear and example files that I downloaded from
% Physionet seem to be inconsistent with the specification. The code below works for
% one example, but your mileage may vary.
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT
% See also BIDS_TSV, BIOPAC_ACQ, BUCN_TXT, EEGSYNTH_TSV, EVENTS_TSV, LIBERTY_CSV, MAUS_TEXTGRID, MOTION_C3D, OPENBCI_TXT, OPENPOSE_KEYPOINTS, OPENSIGNALS_TXT, OPENVIBE_MAT, OPM_FIL, QUALISYS_TSV, SCCN_XDF, SENSYS_CSV, SNIRF, SPIKEGLX_BIN, UNICORN_CSV, XSENS_MVNX

% This code was generated using DeepSeek and subsequently tested and fixed by Robert
% Oostenveld, hence I don't know the appropriate copyrights.
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org for the
% documentation and details.
%
% $Id$

needhdr = (nargin==1);
needevt = (nargin==2);
needdat = (nargin==5);

% Construct file names
[p, f, x] = fileparts(filename);

headerfile = fullfile(p, [f '.hea']);
datafile = fullfile(p, [f '.dat']); % this is what we expect, but it could be different

% Check if files exist
if ~exist(headerfile, 'file')
  error('Header file %s not found', headerfile);
end
if ~exist(datafile, 'file')
  error('Data file %s not found', datafile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(headerfile, 'r');
if fid == -1
  error('Could not open header file %s', headerfile);
end

% Read first line (record line)
first_line = fgetl(fid);
record_info = textscan(first_line, '%s %d %f %d');
num_signals = record_info{2};
Fs = record_info{3};

% Initialize variables
file_name = cell(num_signals, 1);
format = zeros(num_signals, 1);
gain = zeros(num_signals, 1);
baseline = zeros(num_signals, 1);
units = cell(num_signals, 1);
adc_gain = zeros(num_signals, 1);
adc_zero = zeros(num_signals, 1);

% Read signal specifications
for i = 1:num_signals
  signal_line = fgetl(fid);
  if signal_line == -1
    error('Unexpected end of header file');
  end

  % Parse signal line
  parts = textscan(signal_line, '%s %d %s %d %d %f %d %f %s');

  % Store signal information
  file_name{i}  = parts{1}{1};
  format(i)     = parts{2};
  gain(i)       = parts{6};
  baseline(i)   = parts{8};
  units{i}      = parts{9}{1};
  adc_gain(i)   = parts{6};
  adc_zero(i)   = parts{8};
end

% Read remaining lines as comments
comments = {};
while true
  comment_line = fgetl(fid);
  if comment_line == -1
    break;
  end
  comments{end+1} = comment_line;
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if needhdr
  %% convert the header to FieldTrip format

  hdr.Fs = Fs;
  hdr.label = {};
  for i=1:num_signals
    % create channel labels
    hdr.label{i} = sprintf('%d', i);
  end
  hdr.nChans = num_signals;
  hdr.nSamples = Inf;

  % assume continuous data
  hdr.nSamplesPre = 0;
  hdr.nTrials = 0;

  % also store the original header details
  hdr.orig.file_name = file_name;
  hdr.orig.format = format;
  hdr.orig.gain = gain;
  hdr.orig.baseline = baseline;
  hdr.orig.units = units;
  hdr.orig.adc_gain = adc_gain;
  hdr.orig.adc_zero = adc_zero;
  hdr.orig.comments = comments;

  varargout = {hdr};

elseif needdat
  %% parse the data
  datafile = fullfile(p, file_name{i}); % the actual file name as specified in the header

  fid = fopen(datafile, 'r');
  if fid == -1
    error('Could not open data file %s', datafile);
  end

  % Format	Description
  % 8	      First differences stored as signed 8-bit integers.
  % 16	    16-bit two’s complement integers (little-endian).
  % 24	    24-bit two’s complement integers (little-endian).
  % 32	    32-bit two’s complement integers (little-endian).
  % 61	    16-bit two’s complement integers (big-endian).
  % 80	    8-bit offset binary (unsigned 8-bit, subtract 128 to recover).
  % 160	    16-bit offset binary (unsigned 16-bit, subtract 32,768 to recover).
  % 212	    Packed 12-bit two’s complement samples (compact format, common in PhysioBank).
  % 310	    Packed 10-bit two’s complement samples (legacy format).
  % 311	    Alternative packed 10-bit samples (different packing from 310).
  % 508, 516, 524	Signals compressed with FLAC (8, 16, or 24 bits per sample).

  switch format
    case 8
      % not sure, it might also that these are differences and that a CUMSUM is needed
      signal = fread(fid, [num_signals, inf], 'int8'); 
      fclose(fid);
    case 16
      signal = fread(fid, [num_signals, inf], 'int16');
      fclose(fid);
    case 24
      signal = fread(fid, [num_signals, inf], 'int24');
      fclose(fid);
    case 32
      signal = fread(fid, [num_signals, inf], 'int32');
      fclose(fid);
    case 80
      signal = fread(fid, [num_signals, inf], 'uint8')-128;
      fclose(fid);
    case 160
      signal = fread(fid, [num_signals, inf], 'uint16')-32768;
      fclose(fid);
    otherwise
      ft_error('unsupported subformat (%d)', format);
  end

  % Convert digital values to physical units
  for i = 1:num_signals
    signal(i,:) = (signal(i,:) - adc_zero(i)) / adc_gain(i);
  end

  % I cannot make sense of the format when it has multiple channels
  % convert it to a channels-by-samples with one channel
  varargout = {signal(:)'};

elseif needevt
  %% parse the events
  ft_warning('reading of events is not yet implemented');

  % return the events
  varargout = {[]};

end