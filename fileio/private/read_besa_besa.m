
function [data] = read_besa_besa(filename, header, begsample, endsample, chanindx)
%% Reads BESA .besa format files
% See formatting document <a href="matlab:web(http://www.besa.de/downloads/file-formats/)">here</a>
% 
%
% Use as
%   [header] = read_besa_besa(filename);
% where
%    filename        name of the datafile, including the .besa extension
% This returns a header structure with the following elements
%   header.Fs           sampling frequency
%   header.nChans       number of channels
%   header.nSamples     number of samples per trial
%   header.nSamplesPre  number of pre-trigger samples in each trial
%   header.nTrials      number of trials
%   header.label        cell-array with labels of each channel
%   header.orig         detailled EDF header information
%
% Use as
%   [header] = read_besa_besa(filename, [], chanindx);
% where
%    filename        name of the datafile, including the .edf extension
%    chanindx        index of channels to read (optional, default is all)
%                    Note that since 
% This returns a header structure with the following elements
%   header.Fs           sampling frequency
%   header.nChans       number of channels
%   header.nSamples     number of samples per trial
%   header.nSamplesPre  number of pre-trigger samples in each trial
%   header.nTrials      number of trials
%   header.label        cell-array with labels of each channel
%   header.orig         detailled EDF header information
%
% Or use as
%   [dat] = read_besa_besa(filename, header, begsample, endsample, chanindx);
% where
%    filename        name of the datafile, including the .edf extension
%    header          header structure, see above
%    begsample       index of the first sample to read
%    endsample       index of the last sample to read
%    chanindx        index of channels to read (optional, default is all)
% This returns a Nchans X Nsamples data matrix
% 
% 
% 2016 - Kristopher Anderson, Knight Lab, Helen Wills Neuroscience Institute, University of California, Berkeley


% For debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
warning on;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

switch nargin
  case 1
    chanindx=[];
  case 2
    chanindx=[];
  case 3
    chanindx=begsample;
  case 4
    error('ReadBesaMatlab:ErrorInput','Number of input arguments should be 1,2,3, or 5');
end




needhdr = (nargin==1)||(nargin==3);
needevt = (nargin==2); % Not implemented yet  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
needdat = (nargin==5);

if needhdr
  % Return only header data
  data = read_besa_besa_header(filename);
  % Empty chanindx means we want all channels
  if isempty(chanindx)
    chanindx = 1:data.nChans;
  end
  % Return header data for certain channels
  %  Keep original channel info also
  data.orig.chansel = chanindx; % header.orig.chansel is used for reading channels from data
  data.nChans = numel(chanindx);
  data.label = data.label(chanindx);
  data.orig.channel_info.orig_n_channels = data.orig.channel_info.n_channels;
  data.orig.channel_info.n_channels = numel(chanindx);
  data.orig.channel_info.orig_lsbs = data.orig.channel_info.lsbs;
  data.orig.channel_info.lsbs = data.orig.channel_info.lsbs(chanindx);
  data.orig.channel_info.orig_channel_labels = data.orig.channel_info.channel_labels;
  data.orig.channel_info.channel_labels = data.orig.channel_info.channel_labels(chanindx);
  data.orig.channel_info.orig_channel_states = data.orig.channel_info.channel_states;
  data.orig.channel_info.channel_states = data.orig.channel_info.channel_states(chanindx);
  return
end

% Empty chanindx means we want all channels
if isempty(chanindx)
  chanindx = 1:header.nChans;
end

%% Determine channels to pull from data
% Ignore header.orig.chansel and overwrite with chanindx
%   Note that header (input) and data [dat] returned from this function may not match as a result
channels_to_pull = chanindx;
header.orig.chansel = chanindx;
header.nChans = numel(chanindx);
header.label = header.label(chanindx);
header.orig.channel_info.orig_n_channels = header.orig.channel_info.n_channels;
header.orig.channel_info.n_channels = numel(chanindx);
header.orig.channel_info.orig_lsbs = header.orig.channel_info.lsbs;
header.orig.channel_info.lsbs = header.orig.channel_info.lsbs(chanindx);
header.orig.channel_info.orig_channel_labels = header.orig.channel_info.channel_labels;
header.orig.channel_info.channel_labels = header.orig.channel_info.channel_labels(chanindx);
header.orig.channel_info.orig_channel_states = header.orig.channel_info.channel_states;
header.orig.channel_info.channel_states = header.orig.channel_info.channel_states(chanindx);

%% Open file
[fid,msg] = fopen(filename,'r');
assert(fid~=-1,'ReadBesaMatlab:ErrorOpeningFile',msg);

% Get length of file
fseek(fid,0,'eof');
file_length = ftell(fid);
fseek(fid,0,'bof');

%% Data blocks
if needdat
  % Collect data block info
  data_block_offsets   = header.orig.tags.tags.position(strcmp(header.orig.tags.tags.type,'BDAT'));
  data_block_n_samples = header.orig.tags.tags.n_samples(strcmp(header.orig.tags.tags.type,'BDAT'));
  data_block_samples_beg = ones(numel(data_block_n_samples),1);
  data_block_samples_end   = ones(numel(data_block_n_samples),1)*data_block_n_samples(1);
  for block_n = 2:numel(data_block_n_samples)
    data_block_samples_beg(block_n) = data_block_samples_end(block_n-1) + 1;
    data_block_samples_end(block_n)   = data_block_samples_end(block_n-1) + data_block_n_samples(block_n);
  end
  
  % Choose blocks that contain requested samples
  if isempty(begsample) || begsample < 1
    begsample = 1;
  end
  if isempty(endsample) || endsample > sum(data_block_n_samples)
    endsample = sum(data_block_n_samples);
  end
  blocks_to_pull = [];
  for block_n = 1:numel(data_block_n_samples)
    if( ~(data_block_samples_beg(block_n)<begsample && data_block_samples_end(block_n)<begsample) && ...
        ~(data_block_samples_beg(block_n)>endsample && data_block_samples_end(block_n)>endsample) )
      blocks_to_pull(end+1) = block_n; %#ok<AGROW>
    end
  end
  % These values correspond to indices in each block
  samples_to_pull_beg = ones(numel(data_block_n_samples),1);
  samples_to_pull_end = data_block_n_samples';
  for block_n = blocks_to_pull
    if(data_block_samples_beg(block_n)<begsample)
      samples_to_pull_beg(block_n) = begsample-data_block_samples_beg(block_n)+1;
    end
    if(data_block_samples_end(block_n)>endsample)
      samples_to_pull_end(block_n) = endsample-data_block_samples_beg(block_n)+1;
    end
  end
  % These values will correspond to the output sample number
  data_block_samples_beg = max(data_block_samples_beg-begsample+1, 1);
  data_block_samples_end = min(data_block_samples_end-begsample+1, endsample-begsample+1);
  
  % Check for necessary values
  if(~isfield(header.orig.channel_info,'n_channels'))
    fclose(fid);
    error('ReadBesaMatlab:ErrorNoNChannels','header.orig.channel_info.n_channels does not exist. This is needed for reading data blocks');
  end
  if(~isfield(header.orig.channel_info,'lsbs'))
    % No least significant bit values found, so setting them all to 1.0
    header.orig.channel_info.lsbs = ones(header.orig.channel_info.n_channels,1,'double');
  end
  
  % Loop over all data blocks and add to alldata matrix
  data = zeros(numel(channels_to_pull), endsample-begsample+1);
  for block_n = blocks_to_pull
    blockdata = read_BDAT(fid, file_length, data_block_offsets(block_n), header.orig.channel_info.orig_n_channels, header.orig.channel_info.orig_lsbs);
    data(:,data_block_samples_beg(block_n):data_block_samples_end(block_n)) = ...
      blockdata(samples_to_pull_beg(block_n):samples_to_pull_end(block_n),channels_to_pull)';
  end
  % NEED TO IMPLEMENT OVERWRITES %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
end


%% DATA BLOCK FUNCTIONS

function block_data = read_BDAT(fid, file_length, bdat_offset, n_channels, lsbs)
%% Read data block
%
% fif - file ID
% file_length - total length of the file in bytes
% bdat_offset - The location of this data block in the file
% n_channels - number of channels
% lsbs - [1 x n_channels array] - int data is multiplied by this for scaling

if min(lsbs < 0)
  lsbs = ones(size(lsbs));
end
% Please note that zero or negative LSB values are not allowed. If a non-positive value is found in the array, a value of "1.f" is used instead.

% Constants for data type
CONST_FLOAT = 0;
CONST_INT16 = 1;
CONST_COMPRESSED = 1;
CONST_UNCOMPRESSED = 0;
% Constants for prefix byte
CONST_NOCOMPRESSION = 0;
CONST_FIRSTSCHEME = 3;
CONST_SECONDSCHEME = 4;
CONST_THIRDSCHEME = 5;
CONST_NOCOMPRESSION_FIRST2INT = 6;
CONST_FIRSTSCHEME_FIRST2INT = 7;
CONST_NOCOMPRESSION_ALLINT = 8;
CONST_ZLIB_DD = 9;
CONST_ZLIB_FIRSTSCHEME = 13;
CONST_ZLIB_SECONDSCHEME = 14;
CONST_ZLIB_THIRDSCHEME = 15;
CONST_ZLIB_FIRSTSCHEME_FIRST2INT = 17;
CONST_ZLIB_SECONDSCHEME_FIRST2INT = 18;
CONST_ZLIB_THIRDSCHEME_FIRST2INT = 19;
CONST_ZLIB_DD_ALLINT = 29;

% Skip to start of BDAT section
if(fseek(fid,double(bdat_offset),'bof') == -1) % double() because Windows can't seek to uint64
  fclose(fid);
  error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BDAT]',bdat_offset);
end

% Read BDAT tag and offset
[~,ofst_BDAT] = read_tag_offset_pair(fid,'BDAT');

% Check that file is not shorter than expected
if(file_length < (ftell(fid)+ofst_BDAT))
  expected_length = ftell(fid)+ofst_BDAT;
  fclose(fid);
  error('ReadBesaMatlab:ErrorFileTooShortForDataBlock','Data block expected file at least %d bytes long but file is %d bytes long',expected_length,file_length);
end

% Determine type of data in this block
read_tag_offset_pair(fid,'DATT');
flag_BDAT = fread(fid,1,'*uint32');
if bitand(flag_BDAT,uint32(1),'uint32')
  data_type = CONST_INT16;
else
  data_type = CONST_FLOAT;
end
if bitand(flag_BDAT,uint32(hex2dec('0010')),'uint32')
  data_comp = CONST_COMPRESSED;
else
  data_comp = CONST_UNCOMPRESSED;
end

% Determine number of samples in this block
read_tag_offset_pair(fid,'DATS');
n_samples = fread(fid,1,'*uint32');

% Read DATA tag and offset
[~,data_block_length] = read_tag_offset_pair(fid,'DATA');
data_block_offset = double(ftell(fid));

% Read data
block_data = zeros(n_samples,n_channels,'double');
switch num2str([data_type data_comp])
  
  case num2str([CONST_INT16 CONST_UNCOMPRESSED])
    % Read int16s, reshape to [n_samples x n_channels], multiply each channel by LSB
    block_data = bsxfun(@times,lsbs', ...
      double(reshape(fread(fid,n_samples*n_channels,'*int16'),[n_samples,n_channels])));
    
  case num2str([CONST_FLOAT CONST_UNCOMPRESSED])
    % Read singles, reshape to [n_samples x n_channels]
    block_data = double(reshape(fread(fid,n_samples*n_channels,'*single'),[n_samples,n_channels]));
    
  case {num2str([CONST_FLOAT CONST_COMPRESSED]),num2str([CONST_INT16 CONST_COMPRESSED])}
    % Compressed data
    for channel_n = 1:n_channels
      prefix_val = fread(fid,1,'*uint8');
      switch prefix_val
        case CONST_NOCOMPRESSION
          % No zlib. No pre-compression. Yes double difference
          % First two elements are int16. Rest are int16.
          block_data(:,channel_n) = double(fread(fid,n_samples,'*int16'));
          % Integrate twice
          block_data(2:end,channel_n) = cumsum(block_data(2:end,channel_n),1);
          block_data(:,channel_n) = cumsum(block_data(:,channel_n),1);
        case CONST_FIRSTSCHEME
          % No zlib. Yes pre-compression. Yes double difference
          % First two elements are int16. Rest are int8 (pre-compressed, first scheme)
          first_2_vals = fread(fid,2,'*int16');
          block_data(:,channel_n) = double(decode_firstscheme(fid,fread(fid,data_block_length-4,'*uint8'), n_samples, first_2_vals));
        case CONST_SECONDSCHEME
          % No zlib. Yes pre-compression. Yes double difference
          % First two elements are int16. Rest are int8 (pre-compressed, second scheme)
          first_2_vals = fread(fid,2,'*int16');
          block_data(:,channel_n) = double(decode_secondscheme(fid,fread(fid,data_block_length-4,'*uint8'), n_samples, first_2_vals));
        case CONST_THIRDSCHEME
          % No zlib. Yes pre-compression. Yes double difference
          % First two elements are int16. Rest are int8 (pre-compressed, third scheme)
          first_2_vals = fread(fid,2,'*int16');
          block_data(:,channel_n) = double(decode_thirdscheme(fid,fread(fid,data_block_length-4,'*uint8'), n_samples, first_2_vals));
        case CONST_NOCOMPRESSION_FIRST2INT
          % No zlib. No pre-compression. Yes double difference
          % First two elements are int32. Rest are int16
          block_data(1:2,channel_n) = double(fread(fid,2,'*int32'));
          block_data(3:end,channel_n) = double(fread(fid,n_samples-2,'*int16'));
          % Integrate twice
          block_data(2:end,channel_n) = cumsum(block_data(2:end,channel_n),1);
          block_data(:,channel_n) = cumsum(block_data(:,channel_n),1);
        case CONST_FIRSTSCHEME_FIRST2INT
          % No zlib. Yes pre-compression. Yes double difference
          % First two elements are int32. Rest are int8 (pre-compressed, first scheme)
          first_2_vals = fread(fid,2,'*int32');
          block_data(:,channel_n) = double(decode_firstscheme(fid,fread(fid,data_block_length-8,'*uint8'), n_samples, first_2_vals));
        case CONST_NOCOMPRESSION_ALLINT
          % No zlib. No pre-compression. Yes double difference
          % First two elements are int32. Rest are int32
          block_data(:,channel_n) = double(fread(fid,n_samples,'*int32'));
          % Integrate twice
          block_data(2:end,channel_n) = cumsum(block_data(2:end,channel_n),1);
          block_data(:,channel_n) = cumsum(block_data(:,channel_n),1);
        case CONST_ZLIB_DD
          % Yes zlib. No pre-compression. Yes double difference
          % First two elements are int16. Rest are int16
          buffer_len = fread(fid,1,'*uint32');
          buffer_data = fread(fid,buffer_len,'*uint8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          block_data(:,channel_n) = typecast(buffer_data,'int16');
        case CONST_ZLIB_FIRSTSCHEME
          % Yes zlib. Yes pre-compression. Yes double difference
          % First two elements are int16. Rest are int8 (pre-compressed, first scheme)
          buffer_len = fread(fid,1,'*uint32');
          buffer_data = fread(fid,buffer_len,'*uint8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          block_data(:,channel_n) = double(decode_firstscheme(fid,buffer_data(5:end), n_samples, typecast(buffer_data(1:4),'int16')));
        case CONST_ZLIB_SECONDSCHEME
          % Yes zlib. Yes pre-compression. Yes double difference
          % First two elements are int16. Rest are int8 (pre-compressed, second scheme)
          buffer_len = fread(fid,1,'*uint32');
          buffer_data = fread(fid,buffer_len,'*uint8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          block_data(:,channel_n) = double(decode_secondscheme(fid,buffer_data(5:end), n_samples, typecast(buffer_data(1:4),'int16')));
        case CONST_ZLIB_THIRDSCHEME
          % Yes zlib. Yes pre-compression. Yes double difference
          % First two elements are int16. Rest are int8 (pre-compressed, third scheme)
          buffer_len = fread(fid,1,'*uint32');
          buffer_data = fread(fid,buffer_len,'*uint8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          block_data(:,channel_n) = double(decode_thirdscheme(fid,buffer_data(5:end), n_samples, typecast(buffer_data(1:4),'int16')));
        case CONST_ZLIB_FIRSTSCHEME_FIRST2INT
          % Yes zlib. Yes pre-compression. Yes double difference
          % First two elements are int16. Rest are int8 (pre-compressed, first scheme)
          buffer_len = fread(fid,1,'*uint32');
          buffer_data = fread(fid,buffer_len,'*uint8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          block_data(:,channel_n) = double(decode_firstscheme(fid,buffer_data(9:end), n_samples, typecast(buffer_data(1:8),'int32')));
        case CONST_ZLIB_SECONDSCHEME_FIRST2INT
          % Yes zlib. Yes pre-compression. Yes double difference
          % First two elements are int32. Rest are int8 (pre-compressed, second scheme)
          buffer_len = fread(fid,1,'*uint32');
          buffer_data = fread(fid,buffer_len,'*uint8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          block_data(:,channel_n) = double(decode_secondscheme(fid,buffer_data(9:end), n_samples, typecast(buffer_data(1:8),'int32')));
        case CONST_ZLIB_THIRDSCHEME_FIRST2INT
          % Yes zlib. Yes pre-compression. Yes double difference
          % First two elements are int32. Rest are int8 (pre-compressed, third scheme)
          buffer_len = fread(fid,1,'*uint32');
          buffer_data = fread(fid,buffer_len,'*uint8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          block_data(:,channel_n) = double(decode_thirdscheme(fid,buffer_data(9:end), n_samples, typecast(buffer_data(1:8),'int32')));
        case CONST_ZLIB_DD_ALLINT
          % Yes zlib. No pre-compression. Yes double difference
          % First two elements are int32. Rest are int32
          buffer_len = fread(fid,1,'*uint32');
          buffer_data = fread(fid,buffer_len,'*int8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          block_data(:,channel_n) = typecast(buffer_data,'int32'); 
        otherwise
          current_loc = ftell(fid);
          fclose(fid);
          error('ReadBesaMatlab:ErrorBDATReadPrefixValueUnknownScheme','Unknown scheme  CH:%d  prefix_val:%d  File offset:%d',channel_n,prefix_val,current_loc);
      end
    end
    
    if(strcmp(num2str([data_type data_comp]),num2str([CONST_INT16 CONST_COMPRESSED])))
      % Multiply int16 data by lsbs
      block_data = bsxfun(@times,lsbs',block_data);
    end
end

% Check that expected amout of data was read
if((data_block_offset+double(data_block_length)) ~= ftell(fid))
  warning('ReadBesaMatlab:WarningDidNotReadExactBlockLength','%d bytes off. Read %d bytes from data block. Should have read %d bytes', ...
    (ftell(fid)-data_block_offset)-double(data_block_length),ftell(fid)-data_block_offset,double(data_block_length));
end

function outbuffer = decode_firstscheme(fid, inbuffer, n_samples, first2vals)
% Read data in first scheme

CONST_MESHGRID_VALS_1 = -7:7;
CONST_AB_INT32_RANGE = 241:-1:236; % Reverse order. This is needed to determine n_vals
CONST_AB_INT16_RANGE = 247:-1:242;
CONST_AB_INT8_RANGE  = 254:-1:248;

max_lut_val = numel(CONST_MESHGRID_VALS_1)^2-1; % Any buffer value greater than this is an announcing byte

% Use persistent variable so lookup table does not need to be recomputed each time
persistent firstscheme_lookuptable;
if isempty(firstscheme_lookuptable)
  % Create the lookup grid from -7 to 7 in x and y
  [firstscheme_lookuptable(:,:,1),firstscheme_lookuptable(:,:,2)] = meshgrid(CONST_MESHGRID_VALS_1,CONST_MESHGRID_VALS_1);
  % Reshape the lookup grid to be [1:225 x 1:2]
  firstscheme_lookuptable = reshape(firstscheme_lookuptable,[numel(CONST_MESHGRID_VALS_1)^2 2]);
end

% Initialize outbuffer
outbuffer = zeros(n_samples,1,'int32');

% Fill in the first two values
outbuffer(1:2) = first2vals;

% Find first announcing byte (AB) (value outside of LUT)
ab_idx = find(inbuffer>max_lut_val,1,'first');

last_outbuffer_idx = 2; % first2vals
if isempty(ab_idx)
  % No ABs, just use lookup table for whole inbuffer
  % Get the output from the lookup table
  %   Transpose and then use linear indexing in the output to put all
  %   elements into a 1-d array
  try
    outbuffer((last_outbuffer_idx+1):end) = firstscheme_lookuptable(inbuffer+1); % Plus 1 because indices start at 0
  catch ME
    if(strcmp(ME.identifier,'MATLAB:subsassignnumelmismatch'))
      expected_samples = numel(outbuffer((last_outbuffer_idx+1):end));
      received_samples = numel(firstscheme_lookuptable(inbuffer+1));
      fclose(fid);
      error('ReadBesaMatlab:ErrorUnexpectedNSamplesFromPreCompression','Expected %d samples, but got %d samples. [first scheme, no ABs]', ...
        expected_samples,received_samples);
    else
      rethrow(ME);
    end
  end
end

% Loop until out of announcing bytes
possible_abs = inbuffer > max_lut_val;
last_ab_idx = 0;
while ~isempty(ab_idx)
  
  % Fill outbuffer using LUT with all values between the last set of non-encodable values 
  %   and the current set of non-encodable values,
  %   starting at the last filled outbuffer index.
  try
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+2*(ab_idx-last_ab_idx-1))) = ...
      firstscheme_lookuptable(inbuffer((last_ab_idx+1):(ab_idx-1))+1,:); % Plus 1 because indices start at 0
  catch ME
    if(strcmp(ME.identifier,'MATLAB:subsassignnumelmismatch'))
      expected_samples = numel(outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+2*(ab_idx-last_ab_idx-1))));
      received_samples = numel(firstscheme_lookuptable(inbuffer((last_ab_idx+1):(ab_idx-1))+1,:));
      fclose(fid);
      error('ReadBesaMatlab:ErrorUnexpectedNSamplesFromPreCompression','Expected %d samples, but got %d samples. [first scheme, middle of buffer]', ...
        expected_samples,received_samples);
    else
      rethrow(ME);
    end
  end
  last_outbuffer_idx = (last_outbuffer_idx+2*(ab_idx-last_ab_idx-1));
  
  if(any(CONST_AB_INT32_RANGE == inbuffer(ab_idx)))
    % AB indicates int32
    n_vals = find(CONST_AB_INT32_RANGE==inbuffer(ab_idx),1);
    n_skip = n_vals*4; % x4 for int32
    % Fill outbuffer with n_vals
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+n_vals)) = typecast(inbuffer((ab_idx+1):(ab_idx+n_skip)),'int32');
    last_outbuffer_idx = last_outbuffer_idx+n_vals;
    last_ab_idx = ab_idx+n_skip;
  elseif(any(CONST_AB_INT16_RANGE == inbuffer(ab_idx)))
    % AB indicates int16
    n_vals = find(CONST_AB_INT16_RANGE==inbuffer(ab_idx),1);
    n_skip = n_vals*2; % x2 for int16
    % Fill outbuffer with n_vals
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+n_vals)) = typecast(inbuffer((ab_idx+1):(ab_idx+n_skip)),'int16');
    last_outbuffer_idx = last_outbuffer_idx+n_vals;
    last_ab_idx = ab_idx+n_skip;
  elseif(any(CONST_AB_INT8_RANGE == inbuffer(ab_idx)))
    % AB indicates int8
    n_vals = find(CONST_AB_INT8_RANGE==inbuffer(ab_idx),1);
    n_skip = n_vals; % x1 for int8
    % Fill outbuffer with n_vals
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+n_vals)) = typecast(inbuffer((ab_idx+1):(ab_idx+n_skip)),'int8');
    last_outbuffer_idx = last_outbuffer_idx+n_vals;
    last_ab_idx = ab_idx+n_skip;
  else
    % not an alowed announcing byte value
    fclose(fid);
    error('ReadBesaMatlab:ErrorABOutOfRange','Announcing byte out of range: %d',inbuffer(ab_idx));
  end
  
  % Go to next AB
  ab_idx = last_ab_idx + find(possible_abs((last_ab_idx+1):end),1,'first'); % Note: X+[]=[]
  
end

if(last_ab_idx<numel(inbuffer))
  % Fill outbuffer using LUT with all values after the last set of non-encodable values
  %   starting at the last filled outbuffer index.
  try
    outbuffer((last_outbuffer_idx+1):end) = ...
      firstscheme_lookuptable(inbuffer((last_ab_idx+1):end)+1,:); % Plus 1 because indices start at 0
  catch ME
    if(strcmp(ME.identifier,'MATLAB:subsassignnumelmismatch'))
      expected_samples = numel(outbuffer((last_outbuffer_idx+1):end));
      received_samples = numel(firstscheme_lookuptable(inbuffer((last_ab_idx+1):end)+1,:));
      fclose(fid);
      error('ReadBesaMatlab:ErrorUnexpectedNSamplesFromPreCompression','Expected %d samples, but got %d samples. [first scheme, end of buffer]', ...
        expected_samples,received_samples);
    else
      rethrow(ME);
    end
  end
end

% Integrate twice
outbuffer(2:end) = cumsum(outbuffer(2:end));
outbuffer = cumsum(outbuffer);

function outbuffer = decode_secondscheme(fid, inbuffer, n_samples, first2vals)
% Decode second scheme

CONST_MESHGRID_VALS_2A = -2:2;
CONST_MESHGRID_VALS_2B = -5:5;
CONST_AB_INT16_RANGE = 249:-1:246; % Reverse order. This is needed to determine n_vals
CONST_AB_INT8_RANGE  = 254:-1:250;
meshgrid_vals.A = CONST_MESHGRID_VALS_2A;
meshgrid_vals.B = CONST_MESHGRID_VALS_2B;

max_lut_val = numel(CONST_MESHGRID_VALS_2A)^3 + numel(CONST_MESHGRID_VALS_2B)^2 - 1; % Any buffer value greater than this is an announcing byte

% Initialize outbuffer
outbuffer = zeros(n_samples,1,'int32');

% Fill in the first two values
outbuffer(1:2) = first2vals;

% Find first announcing byte (AB) (value outside of LUT)
ab_idx = find(inbuffer>max_lut_val,1,'first');

last_outbuffer_idx = 2; % first2vals
if isempty(ab_idx)
  % No ABs, just use lookup table for whole inbuffer
  % Get the output from the lookup table
  try
    outbuffer((last_outbuffer_idx+1):end) = secondscheme_lookup(inbuffer+1,meshgrid_vals); % Plus 1 because indices start at 0
  catch ME
    if(strcmp(ME.identifier,'MATLAB:subsassignnumelmismatch'))
      expected_samples = numel(outbuffer((last_outbuffer_idx+1):end));
      received_samples = numel(secondscheme_lookup(inbuffer+1,meshgrid_vals));
      fclose(fid);
      error('ReadBesaMatlab:ErrorUnexpectedNSamplesFromPreCompression','Expected %d samples, but got %d samples. [second scheme, no ABs]', ...
        expected_samples,received_samples);
    else
      rethrow(ME);
    end
  end
end

% Loop until out of announcing bytes
possible_abs = inbuffer > max_lut_val;
last_ab_idx = 0;
while ~isempty(ab_idx)
  
  % Fill outbuffer using LUT with all values between the last set of non-encodable values 
  %   and the current set of non-encodable values,
  %   starting at the last filled outbuffer index.
  % No error checking, because we don't know how long it should be
  decoded_buffer = secondscheme_lookup(inbuffer((last_ab_idx+1):(ab_idx-1))+1,meshgrid_vals); % Plus 1 because indices start at 0
  outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+numel(decoded_buffer))) = ...
    decoded_buffer;
  last_outbuffer_idx = (last_outbuffer_idx+numel(decoded_buffer));
  clear decoded_buffer;
  
  if(any(CONST_AB_INT16_RANGE == inbuffer(ab_idx)))
    % AB indicates int16
    n_vals = find(CONST_AB_INT16_RANGE==inbuffer(ab_idx),1);
    n_skip = n_vals*2; % x2 for int16
    % Fill outbuffer with n_vals
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+n_vals)) = typecast(inbuffer((ab_idx+1):(ab_idx+n_skip)),'int16');
    last_outbuffer_idx = last_outbuffer_idx+n_vals;
    last_ab_idx = ab_idx+n_skip;
  elseif(any(CONST_AB_INT8_RANGE == inbuffer(ab_idx)))
    % AB indicates int8
    n_vals = find(CONST_AB_INT8_RANGE==inbuffer(ab_idx),1);
    n_skip = n_vals; % x1 for int8
    % Fill outbuffer with n_vals
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+n_vals)) = typecast(inbuffer((ab_idx+1):(ab_idx+n_skip)),'int8');
    last_outbuffer_idx = last_outbuffer_idx+n_vals;
    last_ab_idx = ab_idx+n_skip;
  else
    % not an allowed announcing byte value
    fclose(fid);
    error('ReadBesaMatlab:ErrorABOutOfRange','Announcing byte out of range [second scheme]: %d',inbuffer(ab_idx));
  end
  
  % Go to next AB
  ab_idx = last_ab_idx + find(possible_abs((last_ab_idx+1):end),1,'first'); % Note: X+[]=[]
  
end

if(last_ab_idx<numel(inbuffer))
  % Fill outbuffer using LUT with all values after the last set of non-encodable values
  %   starting at the last filled outbuffer index.
  try
    outbuffer((last_outbuffer_idx+1):end) = ...
      secondscheme_lookup(inbuffer((last_ab_idx+1):end)+1,meshgrid_vals); % Plus 1 because indices start at 0
  catch ME
    if(strcmp(ME.identifier,'MATLAB:subsassignnumelmismatch'))
      expected_samples = numel(outbuffer((last_outbuffer_idx+1):end));
      received_samples = numel(secondscheme_lookup(inbuffer((last_ab_idx+1):end)+1,meshgrid_vals));
      fclose(fid);
      error('ReadBesaMatlab:ErrorUnexpectedNSamplesFromPreCompression','Expected %d samples, but got %d samples. [second scheme, end of buffer]', ...
        expected_samples,received_samples);
    else
      rethrow(ME);
    end
  end
end

function output = secondscheme_lookup(input, meshgrid_vals)
% Lookup table for second scheme

% Use persistent variable so lookup table does not need to be recomputed each time
persistent secondscheme_lookuptable;
if isempty(secondscheme_lookuptable)
  
  % Create the lookup grid from -2 to 2 in x, y, z
  [secondscheme_lookuptable_a(:,:,:,1),secondscheme_lookuptable_a(:,:,:,2),secondscheme_lookuptable_a(:,:,:,3)] = ...
    meshgrid(meshgrid_vals.A,meshgrid_vals.A,meshgrid_vals.A);
  % Reshape the lookup grid to be [1:125 x 1:3]
  secondscheme_lookuptable_a = reshape(secondscheme_lookuptable_a,[numel(meshgrid_vals.A)^3 3]);
  % Correct order of x,y,z
  secondscheme_lookuptable_a(:,[1 2 3]) = secondscheme_lookuptable_a(:,[3 1 2]);
  
  % Create the lookup grid from -5 to 5 in x and y
  [secondscheme_lookuptable_b(:,:,1),secondscheme_lookuptable_b(:,:,2)] = meshgrid(meshgrid_vals.B,meshgrid_vals.B);
  % Reshape the lookup grid to be [1:121 x 1:2]
  secondscheme_lookuptable_b = reshape(secondscheme_lookuptable_b,[numel(meshgrid_vals.B)^2 2]);
  
  % Put the lookup tables together in a cell array (because of different sized cells)
  secondscheme_lookuptable = num2cell(secondscheme_lookuptable_a,2);
  secondscheme_lookuptable = [secondscheme_lookuptable; num2cell(secondscheme_lookuptable_b,2)];
  
  clear secondscheme_lookuptable_a;
  clear secondscheme_lookuptable_b;
end

output_cell = secondscheme_lookuptable(input);
output = [output_cell{:}];

function outbuffer = decode_thirdscheme(fid, inbuffer, n_samples, first2vals)
% Decode third scheme

CONST_MESHGRID_VALS_3A = -1:1;
CONST_MESHGRID_VALS_3B = -6:6;
CONST_AB_INT16_RANGE = 251:-1:250; % Reverse order. This is needed to determine n_vals
CONST_AB_INT8_RANGE  = 254:-1:252;
meshgrid_vals.A = CONST_MESHGRID_VALS_3A;
meshgrid_vals.B = CONST_MESHGRID_VALS_3B;

max_lut_val = numel(CONST_MESHGRID_VALS_3A)^4 + numel(CONST_MESHGRID_VALS_3B)^2 - 1; % Any buffer value greater than this is an announcing byte

% Initialize outbuffer
outbuffer = zeros(n_samples,1,'int32');

% Fill in the first two values
outbuffer(1:2) = first2vals;

% Find first announcing byte (AB) (value outside of LUT)
ab_idx = find(inbuffer>max_lut_val,1,'first');

last_outbuffer_idx = 2; % first2vals
if isempty(ab_idx)
  % No ABs, just use lookup table for whole inbuffer
  % Get the output from the lookup table
  try
    outbuffer((last_outbuffer_idx+1):end) = thirdscheme_lookup(inbuffer+1,meshgrid_vals); % Plus 1 because indices start at 0
  catch ME
    if(strcmp(ME.identifier,'MATLAB:subsassignnumelmismatch'))
      expected_samples = numel(outbuffer((last_outbuffer_idx+1):end));
      received_samples = numel(thirdscheme_lookup(inbuffer+1,meshgrid_vals));
      fclose(fid);
      error('ReadBesaMatlab:ErrorUnexpectedNSamplesFromPreCompression','Expected %d samples, but got %d samples. [third scheme, no ABs]', ...
        expected_samples,received_samples);
    else
      rethrow(ME);
    end
  end
end

% Loop until out of announcing bytes
possible_abs = inbuffer > max_lut_val;
last_ab_idx = 0;
while ~isempty(ab_idx)
  
  % Fill outbuffer using LUT with all values between the last set of non-encodable values 
  %   and the current set of non-encodable values,
  %   starting at the last filled outbuffer index.
  % No error checking, because we don't know how long it should be
  decoded_buffer = thirdscheme_lookup(inbuffer((last_ab_idx+1):(ab_idx-1))+1,meshgrid_vals); % Plus 1 because indices start at 0
  outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+numel(decoded_buffer))) = ...
    decoded_buffer;
  last_outbuffer_idx = (last_outbuffer_idx+numel(decoded_buffer));
  clear decoded_buffer;
  
  if(any(CONST_AB_INT16_RANGE == inbuffer(ab_idx)))
    % AB indicates int16
    n_vals = find(CONST_AB_INT16_RANGE==inbuffer(ab_idx),1);
    n_skip = n_vals*2; % x2 for int16
    % Fill outbuffer with n_vals
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+n_vals)) = typecast(inbuffer((ab_idx+1):(ab_idx+n_skip)),'int16');
    last_outbuffer_idx = last_outbuffer_idx+n_vals;
    last_ab_idx = ab_idx+n_skip;
  elseif(any(CONST_AB_INT8_RANGE == inbuffer(ab_idx)))
    % AB indicates int8
    n_vals = find(CONST_AB_INT8_RANGE==inbuffer(ab_idx),1);
    n_skip = n_vals; % x1 for int8
    % Fill outbuffer with n_vals
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+n_vals)) = typecast(inbuffer((ab_idx+1):(ab_idx+n_skip)),'int8');
    last_outbuffer_idx = last_outbuffer_idx+n_vals;
    last_ab_idx = ab_idx+n_skip;
  else
    % not an allowed announcing byte value
    fclose(fid);
    error('ReadBesaMatlab:ErrorABOutOfRange','Announcing byte out of range [third scheme]: %d',inbuffer(ab_idx));
  end
  
  % Go to next AB
  ab_idx = last_ab_idx + find(possible_abs((last_ab_idx+1):end),1,'first'); % Note: X+[]=[]
  
end

if(last_ab_idx<numel(inbuffer))
  % Fill outbuffer using LUT with all values after the last set of non-encodable values
  %   starting at the last filled outbuffer index.
  try
    outbuffer((last_outbuffer_idx+1):end) = ...
      thirdscheme_lookup(inbuffer((last_ab_idx+1):end)+1,meshgrid_vals); % Plus 1 because indices start at 0
  catch ME
    if(strcmp(ME.identifier,'MATLAB:subsassignnumelmismatch'))
      expected_samples = numel(outbuffer((last_outbuffer_idx+1):end));
      received_samples = numel(thirdscheme_lookup(inbuffer((last_ab_idx+1):end)+1,meshgrid_vals));
      fclose(fid);
      error('ReadBesaMatlab:ErrorUnexpectedNSamplesFromPreCompression','Expected %d samples, but got %d samples. [third scheme, end of buffer]', ...
        expected_samples,received_samples);
    else
      rethrow(ME);
    end
  end
end

function output = thirdscheme_lookup(input, meshgrid_vals)
% Lookup table for third scheme

% Use persistent variable so lookup table does not need to be recomputed each time
persistent thirdscheme_lookuptable;
if isempty(thirdscheme_lookuptable)
  
  % Create the lookup grid from -1 to 1 in x, y, z, c
  [thirdscheme_lookuptable_a(:,:,:,:,1),thirdscheme_lookuptable_a(:,:,:,:,2),thirdscheme_lookuptable_a(:,:,:,:,3),thirdscheme_lookuptable_a(:,:,:,:,4)] = ...
    ndgrid(meshgrid_vals.A);
  % Reshape the lookup grid to be [1:81 x 1:4]
  thirdscheme_lookuptable_a = reshape(thirdscheme_lookuptable_a,[numel(meshgrid_vals.A)^4 4]);
  % Correct order of x,y,z,c
  thirdscheme_lookuptable_a(:,[1 2 3 4]) = thirdscheme_lookuptable_a(:,[4 3 2 1]);
  
  % Create the lookup grid from -6 to 6 in x and y
  [thirdscheme_lookuptable_b(:,:,1),thirdscheme_lookuptable_b(:,:,2)] = meshgrid(meshgrid_vals.B,meshgrid_vals.B);
  % Reshape the lookup grid to be [1:169 x 1:2]
  thirdscheme_lookuptable_b = reshape(thirdscheme_lookuptable_b,[numel(meshgrid_vals.B)^2 2]);
  
  % Put the lookup tables together in a cell array (because of different sized cells)
  thirdscheme_lookuptable = num2cell(thirdscheme_lookuptable_a,2);
  thirdscheme_lookuptable = [thirdscheme_lookuptable; num2cell(thirdscheme_lookuptable_b,2)];
  
  clear thirdscheme_lookuptable_a;
  clear thirdscheme_lookuptable_b;
end

output_cell = thirdscheme_lookuptable(input);
output = [output_cell{:}];


%% HELPER FUNCTIONS

function M = dunzip(Z)
% DUNZIP - decompress gzipped stream of bytes
% FORMAT M = dzip(Z)
% Z  -  compressed variable to decompress (uint8 vector)
% M  -  decompressed output
% 
% See also DZIP

% Carefully tested, but no warranty; use at your own risk.
% Michael Kleder, Nov 2005
% Modified by Guillaume Flandin, May 2008

import com.mathworks.mlwidgets.io.InterruptibleStreamCopier
a   = java.io.ByteArrayInputStream(Z);
b   = java.util.zip.InflaterInputStream(a);
isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
c   = java.io.ByteArrayOutputStream;
isc.copyStream(b,c);
M   = c.toByteArray;

function [out_tag, out_offset] = read_tag_offset_pair(fid,expected_tag)
% Read 4 bytes and check if they match expected value
out_tag = fread(fid,4,'*char')';
if(nargin>1)
  % Compare tag with expected tag
  if ~strcmp(expected_tag,out_tag)
    curr_offset = ftell(fid);
    fclose(fid);
    error('ReadBesaMatlab:ErrorTagMismatch','Expecting [%s] but read [%s] at offset %d',expected_tag,out_tag,curr_offset);
  end
end
% Read offset value following tag
out_offset = fread(fid,1,'*uint32');






function [header] = read_besa_besa_header(fname)
%% Reads BESA .besa format header information and skips data
% See formatting document <a href="matlab:web(http://www.besa.de/downloads/file-formats/)">here</a>
% 
% [alldata,file_info,channel_info,tags,events] = readbesa(fname)
% 
% inputs:
%  fname [string] - path to .besa file
% 
% outputs:
%  header [structure] - Header information
% 
% 
% 
% 2016 - Kristopher Anderson, Knight Lab, Helen Wills Neuroscience Institute, University of California, Berkeley

% For debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
warning on;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

%% Open file
[fid,msg] = fopen(fname,'r');
assert(fid~=-1,'ReadBesaMatlab:ErrorOpeningFile',msg);

% Get length of file
fseek(fid,0,'eof');
file_length = ftell(fid);
fseek(fid,0,'bof');

%% Header Block
[~,ofst_BCF1] = read_tag_offset_pair(fid,'BCF1');

% Read data in header block
while ~feof(fid) && ftell(fid) < (8+ofst_BCF1) % 8 for header tag ('BCF1') and header offset (uint32)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'VERS'
      % File version
      header.orig.file_info.besa_file_version = read_chars(fid,current_length);
    case 'OFFM'
      % Index of first 'file main info' block (BFMI)
      BFMI_offset = fread(fid,1,'*int64');
    case 'OFTL'
      % Index of first 'tag list' block (BTAG)
      BTAG_offset = fread(fid,1,'*int64');
    case 'OFBI'
      % Index of first 'channel and location' block (BCAL)
      BCAL_offset = fread(fid,1,'*int64');
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if((ftell(fid)+current_length) <= file_length)
        if(fseek(fid,current_length,'cof') == -1)
          fclose(fid);
          error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in header block [BCF1]))',current_length);
        end
      else
        fclose(fid);
        error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
      end
  end
  
end

% Check for necessary header data
if ~exist('BFMI_offset','var')
  fclose(fid);
  error('ReadBesaMatlab:ErrorNoHeaderBFMI','No BFMI block found in header');
end
if ~exist('BTAG_offset','var')
  fclose(fid);
  error('ReadBesaMatlab:ErrorNoHeaderBTAG','No BTAG block found in header');
end
if ~exist('BCAL_offset','var')
  fclose(fid);
  error('ReadBesaMatlab:ErrorNoHeaderBCAL','No BCAL block found in header');
end

%% 'tag list' blocks
header.orig.tags.next_BTAG_ofst = BTAG_offset;
header.orig.tags.offsets = [];
header.orig.tags.n_tags = 0;
% Keep reading until no more BTAG blocks
while header.orig.tags.next_BTAG_ofst > 0
  header.orig.tags = read_BTAG(fid, file_length, header.orig.tags);
end
header.orig.tags = rmfield(header.orig.tags,'next_BTAG_ofst');

% Check that file is not much shorter than expected
%  This does not take into account length of final block but might still be useful
if(file_length <= header.orig.tags.tags.position(end))
  fclose(fid);
  error('ReadBesaMatlab:ErrorFileTooShort','Expected file at least %d bytes long but file is %d bytes long',header.orig.tags.tags(end).position,file_length);
end

%% 'file main info' blocks
header.orig.file_info.next_BFMI_ofst = BFMI_offset;
header.orig.file_info.offsets = [];
% Keep reading until no more BFMI blocks
while header.orig.file_info.next_BFMI_ofst > 0
  header.orig.file_info = read_BFMI(fid, file_length, header.orig.file_info);
end
header.orig.file_info = rmfield(header.orig.file_info,'next_BFMI_ofst');
% NEED TO IMPLEMENT OVERWRITES %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

%% 'channel and location' blocks
header.orig.channel_info.next_BCAL_ofst = BCAL_offset;
header.orig.channel_info.offsets = [];
% Keep reading until no more BCAL blocks
while header.orig.channel_info.next_BCAL_ofst > 0
  header.orig.channel_info = read_BCAL(fid, file_length, header.orig.channel_info);
end
header.orig.channel_info = rmfield(header.orig.channel_info,'next_BCAL_ofst');
% NEED TO IMPLEMENT OVERWRITES %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

if ~isfield(header.orig.channel_info,'n_channels')
  error('ReadBesaMatlab:ErrorNoHeaderNChannels','Missing number of channels in header [BCAL:CHNR]');
end

% Combine info from channel_info.coord_data and channel_info.channel_states to get actual coordinate data
if(isfield(header.orig.channel_info,'channel_states') && isfield(header.orig.channel_info,'coord_data'))
  for channel_n = 1:header.orig.channel_info.n_channels
    %header.orig.channel_info.channel_locations(channel_n) = [];
    header.orig.channel_info.channel_locations(channel_n).x = NaN;
    header.orig.channel_info.channel_locations(channel_n).y = NaN;
    header.orig.channel_info.channel_locations(channel_n).z = NaN;
    header.orig.channel_info.channel_locations(channel_n).xori = NaN; % Orientation
    header.orig.channel_info.channel_locations(channel_n).yori = NaN;
    header.orig.channel_info.channel_locations(channel_n).zori = NaN;
    header.orig.channel_info.channel_locations(channel_n).x2 = NaN; % Second coil
    header.orig.channel_info.channel_locations(channel_n).y2 = NaN;
    header.orig.channel_info.channel_locations(channel_n).z2 = NaN;
    if( header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_SCALPELECTRODE || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_MAGNETOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_AXIAL_GRADIOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_PLANAR_GRADIOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_MEGREFERENCE )
      header.orig.channel_info.channel_locations(channel_n).x = double(header.orig.channel_info.coord_data(channel_n,1));
      header.orig.channel_info.channel_locations(channel_n).y = double(header.orig.channel_info.coord_data(channel_n,2));
      header.orig.channel_info.channel_locations(channel_n).z = double(header.orig.channel_info.coord_data(channel_n,3));
    end
    if( header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_MAGNETOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_AXIAL_GRADIOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_PLANAR_GRADIOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_MEGREFERENCE )
      header.orig.channel_info.channel_locations(channel_n).xori = double(header.orig.channel_info.coord_data(channel_n,7));
      header.orig.channel_info.channel_locations(channel_n).yori = double(header.orig.channel_info.coord_data(channel_n,8));
      header.orig.channel_info.channel_locations(channel_n).zori = double(header.orig.channel_info.coord_data(channel_n,9));
    end
    if( header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_AXIAL_GRADIOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_PLANAR_GRADIOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_MEGREFERENCE )
      header.orig.channel_info.channel_locations(channel_n).x2 = double(header.orig.channel_info.coord_data(channel_n,4));
      header.orig.channel_info.channel_locations(channel_n).y2 = double(header.orig.channel_info.coord_data(channel_n,5));
      header.orig.channel_info.channel_locations(channel_n).z2 = double(header.orig.channel_info.coord_data(channel_n,6));
    end
    if( header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_MEGREFERENCE )
      if( header.orig.channel_info.channel_locations(channel_n).x2==0 && ...
          header.orig.channel_info.channel_locations(channel_n).y2==0 && ...
          header.orig.channel_info.channel_locations(channel_n).z2==0 )
        header.orig.channel_info.channel_locations(channel_n).x2 = NaN;
        header.orig.channel_info.channel_locations(channel_n).y2 = NaN;
        header.orig.channel_info.channel_locations(channel_n).z2 = NaN;
      end
    end
  end
end

%% Events
% Collect event block info
header.orig.events.offsets = header.orig.tags.tags.position(strcmp(header.orig.tags.tags.type,'BEVT'));
header.orig.events.offsets = sort(header.orig.events.offsets, 'ascend'); % Later blocks overwrite matching events
for block_n = 1:numel(header.orig.events.offsets)
  header.orig.events = read_BEVT(fid, file_length, header.orig.events, header.orig.events.offsets(block_n));
end
% NEED TO IMPLEMENT OVERWRITES %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

%% Reorganize header structure
header.nChans = header.orig.channel_info.n_channels;
if isfield(header.orig.file_info,'s_rate')
  header.Fs = header.orig.file_info.s_rate;
else
  warning('ReadBesaMatlab:WarningMissingHeaderInfo','Missing sample rate in header');
  header.Fs = [];
end
if isfield(header.orig.file_info,'n_samples')
  header.nSamples = header.orig.file_info.n_samples;
else
  warning('ReadBesaMatlab:WarningMissingHeaderInfo','Missing number of samples in header');
  header.nSamples = [];
end

header.nSamplesPre = 0; % Continuous data
header.nTrials     = 1;  % Continuous data

%  Channel labels
if isfield(header.orig.channel_info,'channel_labels')
  header.label = header.orig.channel_info.channel_labels;
else
  warning('ReadBesaMatlab:WarningMissingHeaderInfo','Missing channel labels in header.orig. Creating default channel names');
  for channel_n = 1:header.nChans
    header.label{channel_n} = sprintf('chan%03d', channel_n);
  end
end

%  Channel coordinates
if isfield(header.orig.channel_info,'channel_locations')
  for channel_n = 1:header.nChans
    header.elec.label{channel_n} = header.label{channel_n};
    header.elec.pnt(channel_n,1) = header.orig.channel_info.channel_locations(channel_n).x;
    header.elec.pnt(channel_n,2) = header.orig.channel_info.channel_locations(channel_n).y;
    header.elec.pnt(channel_n,3) = header.orig.channel_info.channel_locations(channel_n).z;
  end
end

function tags = read_BTAG(fid, file_length, tags)
%% Read tag block
% tags [structure] - Existing or blank BTAG structure
%                            Blank needs fields:
%                               next_BTAG_ofst - file offset for BTAG to be read
%                               offsets = []
%                               n_tags = 0
% file_length [scalar] - Length of file in bytes
% fid [scalar] - File identifier

% Skip to start of BTAG section
if(fseek(fid,tags.next_BTAG_ofst,'bof') == -1)
  fclose(fid);
  error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BTAG]',tags.next_BTAG_ofst);
end
tags.offsets(end+1) = tags.next_BTAG_ofst;

% Read BTAG tag and offset
[~,tag_block_length] = read_tag_offset_pair(fid,'BTAG');

% Untagged offset to next BTAG section
tags.next_BTAG_ofst = fread(fid,1,'*int64');

% Loop through all tags in data section
while ftell(fid) < (uint64(tags.offsets(end))+uint64(tag_block_length))
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'TAGE'
      % Tag list entry
      tags.n_tags = tags.n_tags+1;
      tags.tags.type{tags.n_tags} = fread(fid,4,'*char')';
      tags.tags.position(tags.n_tags) = fread(fid,1,'*uint64');
      tags.tags.n_samples(tags.n_tags) = double(fread(fid,1,'*uint32'));
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if((ftell(fid)+current_length) <= file_length)
        if(fseek(fid,current_length,'cof') == -1)
          fclose(fid);
          error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BTAG]))',current_length);
        end
      else
        fclose(fid);
        error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
      end
  end
end

% Check that expected amout of file was read
expected_length = double(tag_block_length) + 8; % 8 for tag and offset
if((tags.offsets(end)+expected_length) ~= ftell(fid))
  warning('ReadBesaMatlab:WarningDidNotReadExactBlockLength','%d bytes off. Read %d bytes from tag block. Should have read %d bytes', ...
    (ftell(fid)-tags.offsets(end))-expected_length,ftell(fid)-tags.offsets(end),expected_length);
end

function file_info = read_BFMI(fid, file_length, file_info)
%% Read file main info block
% file_info [structure] - Existing or blank BFMI structure
%                            Blank needs fields:
%                               next_BFMI_ofst - file offset for BFMI to be read
%                               offsets = []
% file_length [scalar] - Length of file in bytes
% fid [scalar] - File identifier

% Skip to start of BFMI section
if(fseek(fid,file_info.next_BFMI_ofst,'bof') == -1)
  fclose(fid);
  error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BFMI]',file_info.next_BFMI_ofst);
end
file_info.offsets(end+1) = file_info.next_BFMI_ofst;

% Read BFMI tag and offset
[~,fileinfo_block_length] = read_tag_offset_pair(fid,'BFMI');

% Untagged offset to next BFMI section
file_info.next_BFMI_ofst = fread(fid,1,'*int64');

% Create staff field if it doesn't exist already. This is necessary because
%   there is no indication of how many staff to expect, so to increment an
%   array, you need an existing array
if(~isfield(file_info,'staff'))
  file_info.staff = [];
end

% Loop through all tags in data section
while ftell(fid) < (uint64(file_info.offsets(end))+uint64(fileinfo_block_length))
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'SAMT'
      % Total number of samples
      file_info.n_samples = double(fread(fid,1,'*int64'));
    case 'SAMP'
      % Number of samples per second
      file_info.s_rate = fread(fid,1,'*double');
    case 'FINN'
      % Name of the institution
      file_info.institution.name = read_chars(fid,current_length);
    case 'FINA'
      % Address of the institution
      fina_end = ftell(fid)+current_length;
      while ~feof(fid) && ftell(fid) < fina_end
        [current_tag,current_length] = read_tag_offset_pair(fid);
        switch current_tag
          case 'ASTR'
            % Street name
            file_info.institution.street_name = read_chars(fid,current_length);
          case 'ASTA'
            % State
            file_info.institution.state = read_chars(fid,current_length);
          case 'ACIT'
            % City
            file_info.institution.city = read_chars(fid,current_length);
          case 'APOS'
            % Post code
            file_info.institution.post_code = read_chars(fid,current_length);
          case 'ACOU'
            % Country
            file_info.institution.country = read_chars(fid,current_length);
          case 'APHO'
            % Phone number
            file_info.institution.phone_number = read_chars(fid,current_length);
          otherwise
            % Unrecognzed tag. Try to skip forward by offset
            warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
            if((ftell(fid)+current_length) <= file_length)
              if(fseek(fid,current_length,'cof') == -1)
                fclose(fid);
                error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BFMI:FINA]))',current_length);
              end
            else
              fclose(fid);
              error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
            end
        end
      end
      clear fina_end;
    case 'FENA'
      % Encryption algorithm
      file_info.encryption = read_chars(fid,current_length);
    case 'FCOA'
      % Compression algorithm
      file_info.compression = read_chars(fid,current_length);
    case 'RECD'
      % Recording start date and time
      file_info.recording_date.start = read_chars(fid,current_length);
    case 'RECE'
      % Recording end date and time
      file_info.recording_date.end = read_chars(fid,current_length);
    case 'RECO'
      % Recording offset to GMT
      file_info.recording_date.gmt_offset = fread(fid,1,'*single');
    case 'RECS'
      % Recording system
      file_info.recording_system.name = read_chars(fid,current_length);
    case 'RIBN'
      % Name of the input box
      file_info.recording_system.info = read_chars(fid,current_length);
    case 'RESW'
      % Name of recording software
      file_info.recording_system.software = read_chars(fid,current_length);
    case 'RATC'
      % Amplifier time constant
      file_info.recording_system.time_constant = fread(fid,1,'*single');
    case 'RSEQ'
      % Sequence number
      file_info.sequence_n = double(fread(fid,1,'*uint32'));
    case 'RSID'
      % Session unique identifier
      file_info.session_id = read_chars(fid,current_length);
    case 'RSNR'
      % Session number
      file_info.sequence_n = double(fread(fid,1,'*int32'));
    case 'RSTC'
      % Study comment
      file_info.comment = read_chars(fid,current_length);
    case 'RSTA'
      % Responsible staff
      % This assumes that, for each staff member, all fields are contiguous
      %   Otherwise, the indices may not line up
      file_info.staff(end+1).name = '';
      file_info.staff(end+1).initials = '';
      file_info.staff(end+1).function = '';
      rsta_end = ftell(fid)+current_length;
      while ~feof(fid) && ftell(fid) < rsta_end
        [current_tag,current_length] = read_tag_offset_pair(fid);
        switch current_tag
          case 'SNAM'
            % Name
            file_info.staff(end).name = read_chars(fid,current_length);
          case 'ASTA'
            % Initials
            file_info.staff(end).initials = read_chars(fid,current_length);
          case 'ACIT'
            % Function
            file_info.staff(end).function = read_chars(fid,current_length);
          otherwise
            % Unrecognzed tag. Try to skip forward by offset
            warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
            if((ftell(fid)+current_length) <= file_length)
              if(fseek(fid,current_length,'cof') == -1)
                fclose(fid);
                error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BFMI:RSTA]))',current_length);
              end
            else
              fclose(fid);
              error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
            end
        end
      end
      clear rsta_end;
    case 'PNAF'
      % Subject first name
      file_info.subject.name.first = read_chars(fid,current_length);
    case 'PNAM'
      % Subject middle name
      file_info.subject.name.middle = read_chars(fid,current_length);
    case 'PATN'
      % Subject last name
      file_info.subject.name.last = read_chars(fid,current_length);
    case 'PNAA'
      % Anonymized subject name
      file_info.subject.anon_name = read_chars(fid,current_length);
    case 'PNAT'
      % Subject title
      file_info.subject.title = read_chars(fid,current_length);
    case 'PATD'
      % Subject date of birth
      file_info.subject.birthdate = read_chars(fid,current_length);
    case 'PDOD'
      % Subject date of death
      file_info.subject.deathdate = read_chars(fid,current_length);
    case 'PAGE'
      % Subject gender
      file_info.subject.gender = read_chars(fid,current_length);
    case 'PAWE'
      % Subject weight
      file_info.subject.weight = fread(fid,1,'*single');
    case 'PAHE'
      % Subject height
      file_info.subject.height = fread(fid,1,'*single');
    case 'PAMS'
      % Subject marital status
      file_info.subject.marital_status = read_chars(fid,current_length);
    case 'PAAD'
      % Subject address
      paad_end = ftell(fid)+current_length;
      while ~feof(fid) && ftell(fid) < paad_end
        [current_tag,current_length] = read_tag_offset_pair(fid);
        switch current_tag
          case 'ASTR'
            % Street name
            file_info.subject.address.street_name = read_chars(fid,current_length);
          case 'ASTA'
            % State
            file_info.subject.address.state = read_chars(fid,current_length);
          case 'ACIT'
            % City
            file_info.subject.address.city = read_chars(fid,current_length);
          case 'APOS'
            % Post code
            file_info.subject.address.post_code = read_chars(fid,current_length);
          case 'ACOU'
            % Country
            file_info.subject.address.country = read_chars(fid,current_length);
          case 'APHO'
            % Phone number
            file_info.subject.address.phone_number = read_chars(fid,current_length);
          otherwise
            % Unrecognzed tag. Try to skip forward by offset
            warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
            if((ftell(fid)+current_length) <= file_length)
              if(fseek(fid,current_length,'cof') == -1)
                fclose(fid);
                error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BFMI:PAAD]))',current_length);
              end
            else
              fclose(fid);
              error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
            end
        end
      end
      clear paad_end;
    case 'PALA'
      % Subject language
      file_info.subject.language = read_chars(fid,current_length);
    case 'PAMH'
      % Subject medical history
      file_info.subject.medical_history = read_chars(fid,current_length);
    case 'PATC'
      % Subject comment
      file_info.subject.comment = read_chars(fid,current_length);
    case 'PATI'
      % Subject ID
      file_info.subject.id = read_chars(fid,current_length);
    case 'INF1'
      % Additional information 1
      file_info.additional_info.inf1 = read_chars(fid,current_length);
    case 'INF2'
      % Additional information 2
      file_info.additional_info.inf2 = read_chars(fid,current_length);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if((ftell(fid)+current_length) <= file_length)
        if(fseek(fid,current_length,'cof') == -1)
          fclose(fid);
          error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BFMI]))',current_length);
        end
      else
        fclose(fid);
        error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
      end
  end
end

% Check that expected amout of file was read
expected_length = double(fileinfo_block_length) + 8; % 8 for tag and offset
if((file_info.offsets(end)+expected_length) ~= ftell(fid))
  warning('ReadBesaMatlab:WarningDidNotReadExactBlockLength','%d bytes off. Read %d bytes from file info block. Should have read %d bytes', ...
    (ftell(fid)-file_info.offsets(end))-expected_length,ftell(fid)-file_info.offsets(end),expected_length);
end

function channel_info = read_BCAL(fid, file_length, channel_info)
%% Read channel info block
% channel_info [structure] - Existing or blank BCAL structure
%                            Blank needs fields:
%                               next_BFMI_ofst - file offset for BCAL to be read
%                               offsets = []
% file_length [scalar] - Length of file in bytes
% fid [scalar] - File identifier

% Skip to start of BCAL section
if(fseek(fid,channel_info.next_BCAL_ofst,'bof') == -1)
  fclose(fid);
  error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BCAL]',channel_info.next_BCAL_ofst);
end
channel_info.offsets(end+1) = channel_info.next_BCAL_ofst;

% Read BCAL tag and offset
[~,channel_block_length] = read_tag_offset_pair(fid,'BCAL');

% Untagged offset to next BCAL section
channel_info.next_BCAL_ofst = fread(fid,1,'*int64');

% Loop through all tags in data section
while ftell(fid) < (uint64(channel_info.offsets(end))+uint64(channel_block_length))
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'CHFL'
      % Channel flag
      channel_info.channel_flags.flag = fread(fid,1,'*uint32');
      channel_info.channel_flags.BSA_ELECTRODE_COORDINATES_FROM_LABELS = logical(bitand(channel_info.channel_flags.flag,uint32(hex2dec('0001')),'uint32'));
      channel_info.channel_flags.BSA_SUPPRESS_SPHERE_TO_ELLIPSOID_TRANSFORMATION = logical(bitand(channel_info.channel_flags.flag,uint32(hex2dec('0002')),'uint32'));
      channel_info.channel_flags.BSA_ELECTRODE_COORDINATES_ON_SPHERE = logical(bitand(channel_info.channel_flags.flag,uint32(hex2dec('0004')),'uint32'));
      channel_info.channel_flags.BSA_ADAPT_SPHERICAL_EEG_TO_MEG_COORDS = logical(bitand(channel_info.channel_flags.flag,uint32(hex2dec('0008')),'uint32'));
      channel_info.channel_flags.BSA_SOURCE_CHANNELS_DERIVED_FROM_MEG = logical(bitand(channel_info.channel_flags.flag,uint32(hex2dec('0010')),'uint32'));
    case 'CHTS'
      % Channel type and states of a channel with the specified index
      channel_n = double(fread(fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      channel_info.channel_states(channel_n).flag = fread(fid,1,'*uint32');
      channel_info.channel_states(channel_n).BSA_CHANTYPE_UNDEFINED = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00000000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_POLYGRAPHIC = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00010000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_TRIGGER = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00020000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_CORTICALGRID = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00040000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_INTRACRANIAL = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00080000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_SCALPELECTRODE = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00100000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_MAGNETOMETER = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00200000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_AXIAL_GRADIOMETER = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00400000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_PLANAR_GRADIOMETER = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('01000000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_MEGREFERENCE = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00800000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_NKC_REFERENCE = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('02000000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_CHANSTATE_BAD = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00000001')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANSTATE_REFERENCE = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00000002')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANSTATE_INTERPOLRECORDED = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00000004')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANSTATE_INVISIBLE = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00001000')),'uint32'));
    case 'CHCO'
      % Channel coordinates in mm
      n_channels = current_length / 4 / 9; % Divide by 4 for *single and by 9 for number of elements per channel
      channel_info.coord_data = zeros(n_channels,9,'single');
      for channel_n = 1:n_channels
        channel_info.coord_data(channel_n,:) = fread(fid,9,'*single');
      end
      % More processing done later to obtain actual coordinates
    case 'CHNR'
      % Total number of channels
      channel_info.n_channels = double(fread(fid,1,'*uint16'));
    case 'CHLA'
      % Channel label of a channel with the specified index
      channel_n = double(fread(fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      channel_info.channel_labels{channel_n} = read_chars(fid,current_length-2); % Subtract 2 from offet for channel_n
    case 'CHET'
      % Electrode thickness
      channel_info.electrode_thickness = fread(fid,1,'*single');
    case 'CHSI'
      % Spline interpolation smoothing constant
      channel_info.spline_smoothing_constant = fread(fid,1,'*single');
    case 'CHLS'
      % Least significant bits of data
      channel_info.lsbs = double(fread(fid,current_length/4,'*single'));
      % Please note that zero or negative LSB values are not allowed. If a non-positive value is found in the array, a value of "1.f" is used instead. %%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'CHSF'
      % Sampling frequency
      channel_info.s_rates = fread(fid,current_length/4,'*double');
    case 'HCMM'
      % Head center in mm
      channel_info.head_center.x = fread(fid,1,'*single');
      channel_info.head_center.y = fread(fid,1,'*single');
      channel_info.head_center.z = fread(fid,1,'*single');
    case 'HRMM'
      % Head radius in mm
      channel_info.head_radius = fread(fid,1,'*single');
    case 'FIDC'
      % Fiducial coordinates in mm
      channel_info.fiducial.nasion.x = fread(fid,1,'*single');
      channel_info.fiducial.nasion.y = fread(fid,1,'*single');
      channel_info.fiducial.nasion.z = fread(fid,1,'*single');
      channel_info.fiducial.left_preauricular.x = fread(fid,1,'*single');
      channel_info.fiducial.left_preauricular.y = fread(fid,1,'*single');
      channel_info.fiducial.left_preauricular.z = fread(fid,1,'*single');
      channel_info.fiducial.right_preauricular.x = fread(fid,1,'*single');
      channel_info.fiducial.right_preauricular.y = fread(fid,1,'*single');
      channel_info.fiducial.right_preauricular.z = fread(fid,1,'*single');
    case 'HSPN'
      % Total number of head surface points
      channel_info.n_addn_surf_pnts = double(fread(fid,1,'*int16'));
    case 'HSPC'
      % Head surface point coordinates
      channel_n = double(fread(fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      channel_info.head_surface_points{channel_n}.x = fread(fid,1,'*single');
      channel_info.head_surface_points{channel_n}.y = fread(fid,1,'*single');
      channel_info.head_surface_points{channel_n}.z = fread(fid,1,'*single');
    case 'HSPD'
      % Head surface point labels
      channel_n = double(fread(fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      channel_info.head_surface_points{channel_n}.label = read_chars(fid,current_length-2); % Subtract 2 from offet for channel_n
    case 'CHCU'
      % Channel units
      channel_n = double(fread(fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      channel_info.channel_units{channel_n} = read_chars(fid,current_length-2); % Subtract 2 from offet for channel_n
    case 'CHFI'
      % Filter information
      offset_chfi = ftell(fid);
      offset_end_chfi = int64(offset_chfi)+int64(current_length);
      channel_n = 0;
      while ftell(fid) < offset_end_chfi
        channel_n=channel_n+1;
        filter_object_offset = ftell(fid);
        filter_object_size = fread(fid,1,'*uint32');
        channel_info.filter_info(channel_n).state.state = fread(fid,1,'*uint32');
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_ACTIVE      = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000001')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_ACTIVE     = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000002')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_NOTCH_ACTIVE          = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000004')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_BAND_PASS_ACTIVE      = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000008')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_PRE_LOWCUTOFF_ACTIVE  = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('01000000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_PRE_SUBTRACT_BASELINE = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('02000000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_TYPE_FORWARD    = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_TYPE_ZERO_PHASE = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000010')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_TYPE_BACKWARD   = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000020')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_SLOPE_06DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_SLOPE_12DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000100')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_SLOPE_24DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000200')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_SLOPE_48DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000300')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_TYPE_ZERO_PHASE = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000300')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_TYPE_FORWARD    = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00001000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_TYPE_BACKWARD   = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00002000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_SLOPE_12DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_SLOPE_06DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00010000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_SLOPE_24DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00020000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_SLOPE_48DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00030000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_NOTCH_REMOVE_2ND_HARMONIC = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00100000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_NOTCH_REMOVE_3RD_HARMONIC = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00200000')),'uint32'));
        channel_info.filter_info(channel_n).low_cutoff  = fread(fid,1,'*single');
        channel_info.filter_info(channel_n).high_cutoff = fread(fid,1,'*single');
        channel_info.filter_info(channel_n).notch_freq  = fread(fid,1,'*single');
        channel_info.filter_info(channel_n).notch_width = fread(fid,1,'*single');
        channel_info.filter_info(channel_n).pass_freq   = fread(fid,1,'*single');
        channel_info.filter_info(channel_n).pass_width  = fread(fid,1,'*single');
        if(ftell(fid) ~= filter_object_offset+filter_object_size)
          warning('ReadBesaMatlab:WarningFilterInfoObject','Did not read expected number bytes in filter object [BCAL:CHFI]. Filter information may be incorrect');
        end
      end
      % Check that expected amout of file was read and move to correct position if not
      if(ftell(fid) ~= offset_end_chfi)
        warning('ReadBesaMatlab:WarningFilterInfoBlock','Did not read expected number of bytes in filter info block [BCAL:CHFI]. Filter information may be incorrect. Skipping to next block');
        fseek(fid,offset_end_chfi+1,'bof');
      end
      % Somewhat complicated structure, no test data %%%%%%%%%%%%%%%%%% TODO
    case 'CHNU'
      % Channel numbers
      channel_info.channel_ns = fread(fid,current_length/4,'*int32');
    case 'CHCM'
      % Channel comments
      channel_n = double(fread(fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      channel_info.channel_comments{channel_n} = read_chars(fid,current_length-2); % Subtract 2 from offet for channel_n
    case 'COMC'
      % BESA CTF component. Internal use only
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BCAL:%s]',current_length,current_tag);
      end
    case 'COMH'
      % BESA head transformation. Internal use only
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BCAL:%s]',current_length,current_tag);
      end
    case 'CHSC'
      % BESA spatial components. Internal use only
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BCAL:%s]',current_length,current_tag);
      end
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if((ftell(fid)+current_length) <= file_length)
        if(fseek(fid,current_length,'cof') == -1)
          fclose(fid);
          error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag [BCAL]))',current_length);
        end
      else
        fclose(fid);
        error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
      end
  end
end

% Check that expected amout of file was read
expected_length = double(channel_block_length) + 8; % 8 for tag and offset
if((channel_info.offsets(end)+expected_length) ~= ftell(fid))
  warning('ReadBesaMatlab:WarningDidNotReadExactBlockLength','%d bytes off. Read %d bytes from channel block. Should have read %d bytes', ...
    (ftell(fid)-channel_info.offsets(end))-expected_length,ftell(fid)-channel_info.offsets(end),expected_length);
end


%% EVENT BLOCK FUNCTIONS

function events = read_BEVT(fid, file_length, events, BEVT_offset)
%% Read event block
% BEVT_offset [scalar] - offset of current event block start
% events [structure] - Existing or blank BEVT structure
%                            Blank needs fields:
%                               offsets - sorted array of location of all event blocks
% file_length [scalar] - Length of file in bytes
% fid [scalar] - File identifier

% Skip to start of BEVT section
if(fseek(fid,BEVT_offset,'bof') == -1)
  fclose(fid);
  error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BEVT]',BEVT_offset);
end

% Read BEVT tag and offset
[~,event_block_length] = read_tag_offset_pair(fid,'BEVT');

% Read LIST tag and offset but don't save anything
read_tag_offset_pair(fid,'LIST');

% Now inside of LIST block
% Read HEAD tag - Assuming that it is first tag in LIST block
[~,head_length] = read_tag_offset_pair(fid,'HEAD');
head_offset = ftell(fid);
% Read data in header block
while ~feof(fid) && ftell(fid) < (head_offset+head_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'EVTS'
      events.n_events = double(fread(fid,1,'*uint32'));
    case 'VERS'
      events.version = double(fread(fid,1,'*uint32'));
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if((ftell(fid)+current_length) <= file_length)
        if(fseek(fid,current_length,'cof') == -1)
          fclose(fid);
          error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:HEAD]))',current_length);
        end
      else
        fclose(fid);
        error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
      end
  end
end

% Read all events as structures and put them into a cell array
for event_n = 1:events.n_events
  [current_tag,current_length] = read_tag_offset_pair(fid);
  events.events{event_n} = read_event_tag(fid,current_tag,current_length);
end

% Check that expected amout of file was read
expected_length = double(event_block_length) + 8; % 8 for tag and offset
if((BEVT_offset+expected_length) ~= ftell(fid))
  warning('ReadBesaMatlab:WarningDidNotReadExactBlockLength','%d bytes off. Read %d bytes from event block. Should have read %d bytes', ...
    (ftell(fid)-BEVT_offset)-expected_length,ftell(fid)-BEVT_offset,expected_length);
end

function event_obj = read_event_tag(fid,event_tag,event_length)
% Read an event into a structure
% Create the event object
event_obj.evttype = event_tag;

switch event_tag
  case 'BASE'
    % Base event tag
    event_obj = read_event_tag_base(fid,ftell(fid),event_length,event_obj);
  case 'COMM'
    % Comment event tag
    event_obj = read_event_tag_comm(fid,ftell(fid),event_length,event_obj);
  case 'MARK'
    % Marker event tag
    event_obj = read_event_tag_mark(fid,ftell(fid),event_length,event_obj);
  case 'GENE'
    % Generic event tag
    event_obj = read_event_tag_gene(fid,ftell(fid),event_length,event_obj);
  case 'SEGM'
    % Segment event tag
    event_obj = read_event_tag_segm(fid,ftell(fid),event_length,event_obj);
  case 'ASGM'
    % Average segment start event tag
    event_obj = read_event_tag_asgm(fid,ftell(fid),event_length,event_obj);
  case 'MPS '
    % Multiple pattern search event tag
    %   used by BESA internally
    if(fseek(fid,event_length,'cof') == -1)
      fclose(fid);
      error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [LIST:MPS]',event_length);
    end
  case 'MPSC'
    % Classified multiple pattern search event tag
    %   used by BESA internally
    if(fseek(fid,event_length,'cof') == -1)
      fclose(fid);
      error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [LIST:MPSC]',event_length);
    end
  case 'PATT'
    % Pattern event tag
    event_obj = read_event_tag_patt(fid,ftell(fid),event_length,event_obj);
  case 'TRIG'
    % Trigger event tag
    event_obj = read_event_tag_trig(fid,ftell(fid),event_length,event_obj);
  case 'PAIR'
    % Paired event tag
    event_obj = read_event_tag_pair(fid,ftell(fid),event_length,event_obj);
  case 'ARTI'
    % Artifact event tag
    event_obj = read_event_tag_arti(fid,ftell(fid),event_length,event_obj);
  case 'EPOC'
    % Epoch event tag
    event_obj = read_event_tag_epoc(fid,ftell(fid),event_length,event_obj);
  case 'IMP '
    % Impedance event tag
    event_obj = read_event_tag_imp(fid,ftell(fid),event_length,event_obj);
  otherwise
    % Unrecognzed tag. Try to skip forward by offset
    warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',event_tag,ftell(fid));
    if(fseek(fid,event_length,'cof') == -1)
      fclose(fid);
      error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [LIST]))',event_length);
    end
end

function event_obj = read_event_tag_base(fid,base_offset,base_length, event_obj)
% Read data in the COMM event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (base_offset+base_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'SAMP'
      % Sample index (zero based)
      event_obj.sample_n = fread(fid,1,'*int64');
    case 'TIME'
      % Event time
      event_obj.time.year = fread(fid,1,'*int16');
      event_obj.time.month = fread(fid,1,'*int16');
      event_obj.time.dayOfWeek = fread(fid,1,'*int16');
      event_obj.time.day = fread(fid,1,'*int16');
      event_obj.time.hour = fread(fid,1,'*int16');
      event_obj.time.minute = fread(fid,1,'*int16');
      event_obj.time.second = fread(fid,1,'*int16');
      event_obj.time.milliseconds = fread(fid,1,'*int16');
      event_obj.time.microseconds = fread(fid,1,'*single');
      event_obj.time.stateFlag = fread(fid,1,'*uint64');
    case 'SIDX'
      % Segment index (zero based)
      event_obj.segment_index = fread(fid,1,'*int32');
    case 'CODE'
      % Event code (zero based)
      % This value is used by events of type Pattern (PATT) and Trigger (TRIG)
      %  to store the pattern number and the trigger code. Note that the number/code minus 1 is stored.
      %  Additionally, events of type Artifact (ARTI) and Epoch (EPOC) use the code value
      %  internally due to historical reasons. Other event types may use the code value to
      %  store additional information.
      event_obj.code = fread(fid,1,'*int32');
    case 'EVID'
      % Internal BESA event ID (zero based)
      event_obj.besa_event_id = fread(fid,1,'*int32');
    case 'STAT'
      % Event state
      event_obj.state.value = fread(fid,1,'*uint32');
      event_obj.state.EVT_STATE_MARKED1 = logical(bitand(event_obj.state.value,uint32(hex2dec('00000010')),'uint32'));
      event_obj.state.EVT_STATE_MARKED2 = logical(bitand(event_obj.state.value,uint32(hex2dec('00000020')),'uint32'));
      event_obj.state.EVT_STATE_MARKED3 = logical(bitand(event_obj.state.value,uint32(hex2dec('00000040')),'uint32'));
      event_obj.state.EVT_STATE_MARKED4 = logical(bitand(event_obj.state.value,uint32(hex2dec('00000080')),'uint32'));
      event_obj.state.EVT_STATE_DELETED = logical(bitand(event_obj.state.value,uint32(hex2dec('01000000')),'uint32'));
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:BASE]))',current_length);
      end
  end
end

function event_obj = read_event_tag_comm(fid,comm_offset,comm_length, event_obj)
% Read data in the COMM event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (comm_offset+comm_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'TEXT'
      % Event text
      event_obj.text = read_chars(fid,current_length);
    case 'BASE'
      event_obj = read_event_tag_base(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:COMM]))',current_length);
      end
  end
end

function event_obj = read_event_tag_mark(fid,mark_offset,mark_length, event_obj)
% Read data in the MARK event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (mark_offset+mark_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'BASE'
      event_obj = read_event_tag_base(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:MARK]))',current_length);
      end
  end
end

function event_obj = read_event_tag_gene(fid,gene_offset,gene_length, event_obj)
% Read data in the GENE event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (gene_offset+gene_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'COMM'
      event_obj = read_event_tag_comm(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:GENE]))',current_length);
      end
  end
end

function event_obj = read_event_tag_segm(fid,segm_offset,segm_length, event_obj)
% Read data in the SEGM event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (segm_offset+segm_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'SBEG'
      % Segment start time
      event_obj.segment_start.year = fread(fid,1,'*int16');
      event_obj.segment_start.month = fread(fid,1,'*int16');
      event_obj.segment_start.dayOfWeek = fread(fid,1,'*int16');
      event_obj.segment_start.day = fread(fid,1,'*int16');
      event_obj.segment_start.hour = fread(fid,1,'*int16');
      event_obj.segment_start.minute = fread(fid,1,'*int16');
      event_obj.segment_start.second = fread(fid,1,'*int16');
      event_obj.segment_start.milliseconds = fread(fid,1,'*int16');
      event_obj.segment_start.microseconds = fread(fid,1,'*single');
      event_obj.segment_start.stateFlag = fread(fid,1,'*uint64');
    case 'DAYT'
      % Day time of segment start in microseconds
      event_obj.segment_start.dayt = fread(fid,1,'*double');
    case 'INTE'
      % Sampling interval in microseconds
      event_obj.segment_start.sampling_interval = fread(fid,1,'*double');
    case 'COMM'
      event_obj = read_event_tag_comm(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:SEGM]))',current_length);
      end
  end
end

function event_obj = read_event_tag_asgm(fid,asgm_offset,asgm_length, event_obj)
% Read data in the ASGM event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (asgm_offset+asgm_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'STIM'
      % Prestimulus baseline interval in microseconds
      event_obj.baseline_interval = fread(fid,1,'*double');
    case 'AVRS'
      % Number of averages
      event_obj.n_averages = fread(fid,1,'*int32');
    case 'COMM'
      event_obj = read_event_tag_comm(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:ASGM]))',current_length);
      end
  end
end

function event_obj = read_event_tag_patt(fid,patt_offset,patt_length, event_obj)
% Read data in the PATT event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (patt_offset+patt_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'BASE'
      event_obj = read_event_tag_base(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:PATT]))',current_length);
      end
  end
end

function event_obj = read_event_tag_trig(fid,trig_offset,trig_length, event_obj)
% Read data in the TRIG event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (trig_offset+trig_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'CODE'
      % Event reaction code
      event_obj.reaction_code = fread(fid,1,'*int32');
    case 'TIME'
      % Event reaction time in seconds
      event_obj.reaction_time = fread(fid,1,'*single');
    case 'COMM'
      event_obj = read_event_tag_comm(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:TRIG]))',current_length);
      end
  end
end

function event_obj = read_event_tag_pair(fid,pair_offset,pair_length, event_obj)
% Read data in the PAIR event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (pair_offset+pair_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'PART'
      % Event information of the partner event.
      %   The data section is used to write the partner event information
      %   as a data element (starting with <eventtype> tag ID of partner event).
      [event_tag,event_length] = read_tag_offset_pair(fid);
      switch event_tag
        case 'BASE'
          event_obj.partner_event = read_event_tag_base(fid,ftell(fid),event_length,event_obj);
        case 'COMM'
          event_obj.partner_event = read_event_tag_comm(fid,ftell(fid),event_length,event_obj);
        case 'MARK'
          event_obj.partner_event = read_event_tag_mark(fid,ftell(fid),event_length,event_obj);
        case 'GENE'
          event_obj.partner_event = read_event_tag_gene(fid,ftell(fid),event_length,event_obj);
        case 'SEGM'
          event_obj.partner_event = read_event_tag_segm(fid,ftell(fid),event_length,event_obj);
        case 'ASGM'
          event_obj.partner_event = read_event_tag_asgm(fid,ftell(fid),event_length,event_obj);
        case 'MPS '
          if(fseek(fid,event_length,'cof') == -1)
            fclose(fid);
            error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [LIST:MPS]',event_length);
          end
        case 'MPSC'
          if(fseek(fid,event_length,'cof') == -1)
            fclose(fid);
            error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [LIST:MPSC]',event_length);
          end
        case 'PATT'
          event_obj.partner_event = read_event_tag_patt(fid,ftell(fid),event_length,event_obj);
        case 'TRIG'
          event_obj.partner_event = read_event_tag_trig(fid,ftell(fid),event_length,event_obj);
        case 'PAIR'
          event_obj.partner_event = read_event_tag_pair(fid,ftell(fid),event_length,event_obj);
        case 'ARTI'
          event_obj.partner_event = read_event_tag_arti(fid,ftell(fid),event_length,event_obj);
        case 'EPOC'
          event_obj.partner_event = read_event_tag_epoc(fid,ftell(fid),event_length,event_obj);
        case 'IMP '
          event_obj.partner_event = read_event_tag_imp(fid,ftell(fid),event_length,event_obj);
        otherwise
          % Unrecognzed tag. Try to skip forward by offset
          warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] in PAIR:PART at offset %d',event_tag,ftell(fid));
          if(fseek(fid,event_length,'cof') == -1)
            fclose(fid);
            error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:PAIR:PART]))',event_length);
          end
      end
    case 'COMM'
      event_obj = read_event_tag_comm(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:PAIR]))',current_length);
      end
  end
end

function event_obj = read_event_tag_arti(fid,arti_offset,arti_length, event_obj)
% Read data in the ARTI event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (arti_offset+arti_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'PAIR'
      event_obj = read_event_tag_pair(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:ARTI]))',current_length);
      end
  end
end

function event_obj = read_event_tag_epoc(fid,epoc_offset,epoc_length, event_obj)
% Read data in the EPOC event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (epoc_offset+epoc_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'PAIR'
      event_obj = read_event_tag_pair(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:EPOC]))',current_length);
      end
  end
end

function event_obj = read_event_tag_imp(fid,imp_offset,imp_length, event_obj)
% Read data in the IMP event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (imp_offset+imp_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'FORM'
      % Indicates the format used to store impedance values in VAL
      %  Set to 0 if the impedance status (valid/invalid) is stored
      %  Set to 1 if impedance values (in kOhm) are stored
      event_obj.impedance.format = fread(fid,1,'*int32');
    case 'NR  '
      % Number of channels for which impedance information is stored.
      %   (Number of elements stored in TYPE, LABL and VAL)
      event_obj.impedance.n_channels = fread(fid,1,'*uint32');
    case 'TYPE'
      % Channel types
      % The flags used for channel type description are the same as used for the CHTS data elements (specified in chapter 2.3).
      %   Note: Channel type and channel label are used for identification of channel
      %   for which impedance information is set. (Compare to data elements CHTS and CHLA, specified in chapter 2.3).
      event_obj.impedance.types = fread(fid,event_obj.impedance.n_channels,'*uint32'); % Assumes that n_channels has been set
    case 'LABL'
      % Channel labels
      %   Note: Channel type and channel label are used for identification of channel for which impedance information is set.
      %   (Compare to data elements CHTS and CHLA as specified in chapter 2.3).
      event_obj.impedance.labels = read_chars(fid,current_length);
    case 'VAL '
      % Impedance values
      % Depending on format set in FORM either an impedance STATUS (ok/not ok) or an impedance VALUE (in kOhm) is stored
      %   A value of -1 means that the impedance is not set or invalid
      event_obj.impedance.values = fread(fid,event_obj.impedance.n_channels,'*single'); % Assumes that n_channels has been set
    case 'BASE'
      event_obj = read_event_tag_base(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:IMP]))',current_length);
      end
  end
end

function outchars = read_chars(fid,n_chars)
% Read n_chars characters from file at current position
%   Replace null character (aka char(0) or '\x0') with ''
%   Note transpose after fread
outchars = regexprep(fread(fid,n_chars,'*char')','\x0','');
