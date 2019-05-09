function [varargout] = read_plexon_plx(filename, varargin)

% READ_PLEXON_PLX reads header or data from a Plexon *.plx file, which
% is a file containing action-potential (spike) timestamps and waveforms
% (spike channels), event timestamps (event channels), and continuous
% variable data (continuous A/D channels).
%
% Use as
%   [hdr] = read_plexon_plx(filename)
%   [dat] = read_plexon_plx(filename, ...)
%   [dat1, dat2, dat3, hdr] = read_plexon_plx(filename, ...)
%
% Optional input arguments should be specified in key-value pairs
%   'header'           = structure with header information
%   'memmap'           = 0 or 1
%   'feedback'         = 0 or 1
%   'ChannelIndex'     = number, or list of numbers (that will result in multiple outputs)
%   'SlowChannelIndex' = number, or list of numbers (that will result in multiple outputs)
%   'EventIndex'       = number, or list of numbers (that will result in multiple outputs)

% Copyright (C) 2007-2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% parse the optional input arguments
hdr              = ft_getopt(varargin, 'header');
memmap           = ft_getopt(varargin, 'memmap', false);
feedback         = ft_getopt(varargin, 'feedback', true);
ChannelIndex     = ft_getopt(varargin, 'ChannelIndex');     % type 1
EventIndex       = ft_getopt(varargin, 'EventIndex');       % type 4
SlowChannelIndex = ft_getopt(varargin, 'SlowChannelIndex'); % type 5

needhdr = isempty(hdr);

% start with empty return values
varargout = {};

% the datafile is little endian, hence it may be neccessary to swap bytes in
% the memory mapped data stream depending on the CPU type of this computer
if littleendian
  swapFcn = @(x) x;
else
  swapFcn = @(x) swapbytes(x);
end

% read header info from file, use Matlabs for automatic byte-ordering
fid = fopen_or_error(filename, 'r', 'ieee-le');
fseek(fid, 0, 'eof');
siz = ftell(fid);
fseek(fid, 0, 'bof');

if needhdr
  if feedback, fprintf('reading header from %s\n', filename); end
  % a PLX file consists of a file header, channel headers, and data blocks
  hdr = PL_FileHeader(fid);
  for i=1:hdr.NumDSPChannels
    hdr.ChannelHeader(i) = PL_ChannelHeader(fid);
  end
  for i=1:hdr.NumEventChannels
    hdr.EventHeader(i) = PL_EventHeader(fid);
  end
  for i=1:hdr.NumSlowChannels
    hdr.SlowChannelHeader(i) = PL_SlowChannelHeader(fid);
  end
  hdr.DataOffset = ftell(fid);
  
  if memmap
    % open the file as meory mapped object, note that byte swapping may be needed
    mm = memmapfile(filename, 'offset', hdr.DataOffset, 'format', 'int16');
  end
  
  dum = struct(...
    'Type', [],...
    'UpperByteOf5ByteTimestamp', [],...
    'TimeStamp', [],...
    'Channel', [],...
    'Unit', [],...
    'NumberOfWaveforms', [],...
    'NumberOfWordsInWaveform', [] ...
    );
  
  % read the header of each data block and remember its data offset in bytes
  Nblocks = 0;
  offset  = hdr.DataOffset;  % only used when reading from memmapped file
  hdr.DataBlockOffset = [];
  hdr.DataBlockHeader = dum;
  while offset<siz
    if Nblocks>=length(hdr.DataBlockOffset);
      % allocate another 1000 elements, this prevents continuous reallocation
      hdr.DataBlockOffset(Nblocks+10000) = 0;
      hdr.DataBlockHeader(Nblocks+10000) = dum;
      if feedback, fprintf('reading DataBlockHeader %4.1f%%\n', 100*(offset-hdr.DataOffset)/(siz-hdr.DataOffset)); end
    end
    Nblocks = Nblocks+1;
    if memmap
      % get the header information from the memory mapped file
      hdr.DataBlockOffset(Nblocks) = offset;
      hdr.DataBlockHeader(Nblocks) = PL_DataBlockHeader(mm, offset-hdr.DataOffset, swapFcn);
      % skip the header (16 bytes) and the data (int16 words)
      offset = offset + 16 + 2 * double(hdr.DataBlockHeader(Nblocks).NumberOfWordsInWaveform * hdr.DataBlockHeader(Nblocks).NumberOfWaveforms);
    else
      % read the header information from the file the traditional way
      hdr.DataBlockOffset(Nblocks) = offset;
      hdr.DataBlockHeader(Nblocks) = PL_DataBlockHeader(fid, [], swapFcn);
      fseek(fid, 2 * double(hdr.DataBlockHeader(Nblocks).NumberOfWordsInWaveform * hdr.DataBlockHeader(Nblocks).NumberOfWaveforms), 'cof');  % data consists of short integers
      offset = ftell(fid);
    end % if memmap
  end
  % this prints the final 100%
  if feedback, fprintf('reading DataBlockHeader %4.1f%%\n', 100*(offset-hdr.DataOffset)/(siz-hdr.DataOffset)); end
  % remove the allocated space that was not needed
  hdr.DataBlockOffset = hdr.DataBlockOffset(1:Nblocks);
  hdr.DataBlockHeader = hdr.DataBlockHeader(1:Nblocks);
  
end % if needhdr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the spike channel data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(ChannelIndex)
  if feedback, fprintf('reading spike data from %s\n', filename); end
  
  if memmap
    % open the file as meory mapped object, note that byte swapping may be needed
    mm = memmapfile(filename, 'offset', hdr.DataOffset, 'format', 'int16');
  end
  
  type = [hdr.DataBlockHeader.Type];
  chan = [hdr.DataBlockHeader.Channel];
  ts   = [hdr.DataBlockHeader.TimeStamp];
  
  for i=1:length(ChannelIndex)
    % determine the data blocks with continuous data belonging to this channel
    sel = (type==1 & chan==hdr.ChannelHeader(ChannelIndex(i)).Channel);
    sel = find(sel);
    
    if isempty(sel)
      ft_warning('spike channel %d contains no data', ChannelIndex(i));
      varargin{end+1} = [];
      continue;
    end
    
    % the number of samples can potentially be different in each block
    num    = double([hdr.DataBlockHeader(sel).NumberOfWordsInWaveform]) .* double([hdr.DataBlockHeader(sel).NumberOfWaveforms]);
    
    % check whether the number of samples per block makes sense
    if any(num~=num(1))
      ft_error('spike channel blocks with diffent number of samples');
    end
    
    % allocate memory to hold the data
    buf = zeros(num(1), length(sel), 'int16');
    
    if memmap
      % get the header information from the memory mapped file
      datbeg = double(hdr.DataBlockOffset(sel) - hdr.DataOffset)/2 + 8 + 1;  % expressed in 2-byte words, minus the file header, skip the 16 byte block header
      datend = datbeg + num - 1;
      for j=1:length(sel)
        buf(:,j) = mm.Data(datbeg(j):datend(j));
      end
      % optionally swap the bytes to correct for the endianness
      buf = swapFcn(buf);
    else
      % read the data from the file in the traditional way
      offset = double(hdr.DataBlockOffset(sel)) + 16;  % expressed in bytes, skip the 16 byte block header
      for j=1:length(sel)
        fseek(fid, offset(j), 'bof');
        buf(:,j) = fread(fid, num(j), 'int16');
      end
    end % if memmap
    
    % remember the data for this channel
    varargout{i} = buf;
    
  end %for ChannelIndex
end % if ChannelIndex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the continuous channel data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(SlowChannelIndex)
  if feedback, fprintf('reading continuous data from %s\n', filename); end
  
  if memmap
    % open the file as meory mapped object, note that byte swapping may be needed
    mm = memmapfile(filename, 'offset', hdr.DataOffset, 'format', 'int16');
  end
  
  type = [hdr.DataBlockHeader.Type];
  chan = [hdr.DataBlockHeader.Channel];
  ts   = [hdr.DataBlockHeader.TimeStamp];
  
  for i=1:length(SlowChannelIndex)
    % determine the data blocks with continuous data belonging to this channel
    sel = (type==5 & chan==hdr.SlowChannelHeader(SlowChannelIndex(i)).Channel);
    sel = find(sel);
    
    if isempty(sel)
      ft_error(sprintf('Continuous channel %d contains no data', SlowChannelIndex(i)));
      % ft_warning('Continuous channel %d contains no data', SlowChannelIndex(i));
      % varargin{end+1} = [];
      % continue;
    end
    
    % the number of samples can be different in each block
    num    = double([hdr.DataBlockHeader(sel).NumberOfWordsInWaveform]) .* double([hdr.DataBlockHeader(sel).NumberOfWaveforms]);
    cumnum = cumsum([0 num]);
    
    % allocate memory to hold the data
    buf = zeros(1, cumnum(end), 'int16');
    if memmap
      % get the header information from the memory mapped file
      datbeg = double(hdr.DataBlockOffset(sel) - hdr.DataOffset)/2 + 8 + 1;  % expressed in 2-byte words, minus the file header, skip the 16 byte block header
      datend = datbeg + num - 1;
      for j=1:length(sel)
        bufbeg = cumnum(j)+1;
        bufend = cumnum(j+1);
        % copy the data from the memory mapped file into the continuous buffer
        buf(bufbeg:bufend) = mm.Data(datbeg(j):datend(j));
      end
      % optionally swap the bytes to correct for the endianness
      buf = swapFcn(buf);
    else
      % read the data from the file in the traditional way
      offset = double(hdr.DataBlockOffset(sel)) + 16;  % expressed in bytes, skip the 16 byte block header
      for j=1:length(sel)
        bufbeg = cumnum(j)+1;
        bufend = cumnum(j+1);
        % copy the data from the file into the continuous buffer
        fseek(fid, offset(j), 'bof');
        buf(bufbeg:bufend) = fread(fid, num(j), 'int16');
      end
    end % if memmap
    
    % remember the data for this channel
    varargout{i} = buf;
    
  end %for SlowChannelIndex
end % if SlowChannelIndex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the event channel data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(EventIndex)
  if feedback, fprintf('reading events from %s\n', filename); end
  type = [hdr.DataBlockHeader.Type];
  unit = [hdr.DataBlockHeader.Unit];
  chan = [hdr.DataBlockHeader.Channel];
  ts   = [hdr.DataBlockHeader.TimeStamp];
  
  for i=1:length(EventIndex)
    % determine the data blocks with continuous data belonging to this channel
    sel = (type==4 & chan==hdr.EventHeader(EventIndex(i)).Channel);
    sel = find(sel);
    
    % all information is already contained in the DataBlockHeader, i.e. there is nothing to read
    if isempty(sel)
      ft_warning('event channel %d contains no data', EventIndex(i));
    end
    event.TimeStamp = ts(sel);
    event.Channel   = chan(sel);
    event.Unit      = unit(sel);
    varargout{i}    = event;
    
  end % for EventIndex
end % if EventIndex

fclose(fid);

% always return the header as last
varargout{end+1} = hdr;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS for reading the different header elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = PL_FileHeader(fid)
hdr.MagicNumber         = fread(fid, 1,       'uint32=>uint32');   % = 0x58454c50;
hdr.Version             = fread(fid, 1,       'int32' );           % Version of the data format; determines which data items are valid
hdr.Comment             = fread(fid, [1 128], 'uint8=>char'  );     % User-supplied comment
hdr.ADFrequency         = fread(fid, 1,       'int32' );           % Timestamp frequency in hertz
hdr.NumDSPChannels      = fread(fid, 1,       'int32' );           % Number of DSP channel headers in the file
hdr.NumEventChannels    = fread(fid, 1,       'int32' );           % Number of Event channel headers in the file
hdr.NumSlowChannels     = fread(fid, 1,       'int32' );           % Number of A/D channel headers in the file
hdr.NumPointsWave       = fread(fid, 1,       'int32' );           % Number of data points in waveform
hdr.NumPointsPreThr     = fread(fid, 1,       'int32' );           % Number of data points before crossing the threshold
hdr.Year                = fread(fid, 1,       'int32' );           % Time/date when the data was acquired
hdr.Month               = fread(fid, 1,       'int32' );
hdr.Day                 = fread(fid, 1,       'int32' );
hdr.Hour                = fread(fid, 1,       'int32' );
hdr.Minute              = fread(fid, 1,       'int32' );
hdr.Second              = fread(fid, 1,       'int32' );
hdr.FastRead            = fread(fid, 1,       'int32' );           % reserved
hdr.WaveformFreq        = fread(fid, 1,       'int32' );           % waveform sampling rate; ADFrequency above is timestamp freq
hdr.LastTimestamp       = fread(fid, 1,       'double');           % duration of the experimental session, in ticks
% The following 6 items are only valid if Version >= 103
hdr.Trodalness          = fread(fid, 1,       'char'  );           % 1 for single, 2 for stereotrode, 4 for tetrode
hdr.DataTrodalness      = fread(fid, 1,       'char'  );           % trodalness of the data representation
hdr.BitsPerSpikeSample  = fread(fid, 1,       'char'  );           % ADC resolution for spike waveforms in bits (usually 12)
hdr.BitsPerSlowSample   = fread(fid, 1,       'char'  );           % ADC resolution for slow-channel data in bits (usually 12)
hdr.SpikeMaxMagnitudeMV = fread(fid, 1,       'uint16');           % the zero-to-peak voltage in mV for spike waveform adc values (usually 3000)
hdr.SlowMaxMagnitudeMV  = fread(fid, 1,       'uint16');           % the zero-to-peak voltage in mV for slow-channel waveform adc values (usually 5000); Only valid if Version >= 105 (usually either 1000 or 500)
% The following item is only valid if Version >= 105
hdr.SpikePreAmpGain     = fread(fid, 1,       'uint16');           % so that this part of the header is 256 bytes
hdr.Padding             = fread(fid, 46,      'char'  );           % so that this part of the header is 256 bytes
% Counters for the number of timestamps and waveforms in each channel and unit.
% Note that these only record the counts for the first 4 units in each channel.
% channel numbers are 1-based - array entry at [0] is unused
hdr.TSCounts            = fread(fid, [5 130], 'int32' );            % number of timestamps[channel][unit]
hdr.WFCounts            = fread(fid, [5 130], 'int32' );            % number of waveforms[channel][unit]
%  Starting at index 300, the next array also records the number of samples for the
%  continuous channels.  Note that since EVCounts has only 512 entries, continuous
%  channels above channel 211 do not have sample counts.
hdr.EVCounts            = fread(fid, 512,     'int32' );            % number of timestamps[event_number]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = PL_ChannelHeader(fid)
hdr.Name             = fread(fid, [1 32],  'uint8=>char'   );  % Name given to the DSP channel
hdr.SIGName          = fread(fid, [1 32],  'uint8=>char'   );  % Name given to the corresponding SIG channel
hdr.Channel          = fread(fid, 1,       'int32'  );        % DSP channel number, 1-based
hdr.WFRate           = fread(fid, 1,       'int32'  );        % When MAP is doing waveform rate limiting, this is limit w/f per sec divided by 10
hdr.SIG              = fread(fid, 1,       'int32'  );        % SIG channel associated with this DSP channel 1 - based
hdr.Ref              = fread(fid, 1,       'int32'  );        % SIG channel used as a Reference signal, 1- based
hdr.Gain             = fread(fid, 1,       'int32'  );        % actual gain divided by SpikePreAmpGain. For pre version 105, actual gain divided by 1000.
hdr.Filter           = fread(fid, 1,       'int32'  );        % 0 or 1
hdr.Threshold        = fread(fid, 1,       'int32'  );        % Threshold for spike detection in a/d values
hdr.Method           = fread(fid, 1,       'int32'  );        % Method used for sorting units, 1 - boxes, 2 - templates
hdr.NUnits           = fread(fid, 1,       'int32'  );        % number of sorted units
hdr.Template         = fread(fid, [64 5],  'int16'  );        % Templates used for template sorting, in a/d values
hdr.Fit              = fread(fid, 5,       'int32'  );        % Template fit
hdr.SortWidth        = fread(fid, 1,       'int32'  );        % how many points to use in template sorting (template only)
hdr.Boxes            = reshape(fread(fid, 4*2*5,   'int16'  ), [4 2 5]);        % the boxes used in boxes sorting
hdr.SortBeg          = fread(fid, 1,       'int32'  );        % beginning of the sorting window to use in template sorting (width defined by SortWidth)
hdr.Comment          = fread(fid, [1 128], 'uint8=>char'   );
hdr.Padding          = fread(fid, 11,      'int32'  );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = PL_EventHeader(fid)
hdr.Name          = fread(fid, [1 32],   'uint8=>char'  ); % name given to this event
hdr.Channel       = fread(fid, 1,        'int32'       );       % event number, 1-based
hdr.Comment       = fread(fid, [1 128],  'uint8=>char'  );
hdr.Padding       = fread(fid, 33,       'int32'       );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = PL_SlowChannelHeader(fid)
hdr.Name          = fread(fid, [1 32],  'uint8=>char'  ); % name given to this channel
hdr.Channel       = fread(fid, 1,       'int32' ); % channel number, 0-based
hdr.ADFreq        = fread(fid, 1,       'int32' ); % digitization frequency
hdr.Gain          = fread(fid, 1,       'int32' ); % gain at the adc card
hdr.Enabled       = fread(fid, 1,       'int32' ); % whether this channel is enabled for taking data, 0 or 1
hdr.PreAmpGain    = fread(fid, 1,       'int32' ); % gain at the preamp
% As of Version 104, this indicates the spike channel (PL_ChannelHeader.Channel) of
% a spike channel corresponding to this continuous data channel.
% <=0 means no associated spike channel.
hdr.SpikeChannel  = fread(fid, 1,       'int32' );
hdr.Comment       = fread(fid, [1 128], 'uint8=>char'  );
hdr.Padding       = fread(fid, 28,      'int32' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = PL_DataBlockHeader(fid, offset, swapFcn)
%   % this is the conventional code, it has been replaced by code that works
%   % with both regular and memmapped files
%   hdr.Type                        = fread(fid, 1,  'int16=>int16'   ); % Data type; 1=spike, 4=Event, 5=continuous
%   hdr.UpperByteOf5ByteTimestamp   = fread(fid, 1,  'uint16=>uint16' ); % Upper 8 bits of the 40 bit timestamp
%   hdr.TimeStamp                   = fread(fid, 1,  'uint32=>uint32' ); % Lower 32 bits of the 40 bit timestamp
%   hdr.Channel                     = fread(fid, 1,  'int16=>int16'   ); % Channel number
%   hdr.Unit                        = fread(fid, 1,  'int16=>int16'   ); % Sorted unit number; 0=unsorted
%   hdr.NumberOfWaveforms           = fread(fid, 1,  'int16=>int16'   ); % Number of waveforms in the data to folow, usually 0 or 1
%   hdr.NumberOfWordsInWaveform     = fread(fid, 1,  'int16=>int16'   ); % Number of samples per waveform in the data to follow
if isa(fid, 'memmapfile')
  mm = fid;
  datbeg = offset/2 + 1; % the offset is in bytes (minus the file header), the memory mapped file is indexed in int16 words
  datend = offset/2 + 8;
  buf = mm.Data(datbeg:datend);
else
  buf = fread(fid, 8, 'int16=>int16');
end
hdr.Type                        = swapFcn(buf(1));
hdr.UpperByteOf5ByteTimestamp   = swapFcn(uint16(buf(2)));
hdr.TimeStamp                   = swapFcn(typecast(buf([3 4]), 'uint32'));
hdr.Channel                     = swapFcn(buf(5));
hdr.Unit                        = swapFcn(buf(6));
hdr.NumberOfWaveforms           = swapFcn(buf(7));
hdr.NumberOfWordsInWaveform     = swapFcn(buf(8));
