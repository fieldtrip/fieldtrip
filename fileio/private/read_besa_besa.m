
function [header,data] = read_besa_besa(fname)
%% Reads BESA .besa format files
% See formatting document <a href="matlab:web(http://www.besa.de/downloads/file-formats/)">here</a>
% 
% [header,data] = readbesa(fname)
% 
% inputs:
%  fname [string] - path to .besa file
% 
% outputs:
%  alldata [n_samples x n_channels] - The data
%  header [structure] - Header structure
% 
% 
% 
% 2015 - Kristopher Anderson, Knight Lab, Helen Wills Neuroscience Institute, University of California, Berkeley

% For debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
warning on;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

%% Open file
[fid,msg] = fopen(fname,'r');
assert(fid~=-1,'ReadBesaMatlab:ErrorOpeningFile',msg);

% Get length of file
fseek(fid,0,'eof');
file_length = ftell(fid);
fseek(fid,0,'bof');

%% Read and save all header information
[header] = read_besa_besa_header(fname);

%% Data blocks
% Collect data block info
data_block_offsets   = header.orig.tags.tags.position(strcmp(header.orig.tags.tags.type,'BDAT'));
data_block_n_samples = header.orig.tags.tags.n_samples(strcmp(header.orig.tags.tags.type,'BDAT'));
data_block_samples_begin = ones(numel(data_block_n_samples),1);
data_block_samples_end   = ones(numel(data_block_n_samples),1)*data_block_n_samples(1);
for block_n = 2:numel(data_block_n_samples)
  data_block_samples_begin(block_n) = data_block_samples_end(block_n-1) + 1;
  data_block_samples_end(block_n)   = data_block_samples_end(block_n-1) + data_block_n_samples(block_n);
end

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
data = zeros(sum(data_block_n_samples),header.orig.channel_info.n_channels);
for block_n = 1:numel(data_block_n_samples)
  data(data_block_samples_begin(block_n):data_block_samples_end(block_n),:) = ...
    read_BDAT(fid, file_length, data_block_offsets(block_n), header.orig.channel_info.n_channels, header.orig.channel_info.lsbs);
end
% NEED TO IMPLEMENT OVERWRITES %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO



%% DATA BLOCK FUNCTIONS

function block_data = read_BDAT(fid, file_length, bdat_offset, n_channels, lsbs)
%% Read data block
%
% lsbs - [1 x n_channels array] - int data is multiplied by this for scaling

% CHECK FOR LSB LESS THAN ZERO %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
% Please note that zero or negative LSB values are not allowed. If a non-positive value is found in the array, a value of "1.f" is used instead. %%%%%%%%%%%%%%%%%%%%%%%%%%%

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
      % CHECK FOR LSB LESS THAN ZERO %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
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













