function [dat] = read_neuralynx_dma(filename, begsample, endsample, channel)

% READ_NEURALYNX_DMA reads specified samples and channels data from a Neuralynx DMA log file
%
% Use as
%    [hdr] = read_neuralynx_dma(filename)
%    [dat] = read_neuralynx_dma(filename, begsample, endsample)
%    [dat] = read_neuralynx_dma(filename, begsample, endsample, chanindx)
%
% The channel specification can be a vector with indices, or a single string with the value
%    'all', 'stx', 'pid', 'siz', 'tsh', 'tsl',
%    'cpu', 'ttl', 'x01',  ...,  'x10'
%
% This function returns the electrophysiological data in AD units
% and not in uV. You should look up the details of the headstage and
% the Neuralynx amplifier and scale the values accordingly.

% Copyright (C) 2005-2008, Robert Oostenveld
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The data is simply a stream of constant size records,  the size of the
% record is dependant on the number of Input Boards in the system.  8 in your
% case, there are 32 - 32 bit integers from each Input Board.
%
% // Digital Lynx Packets look like this for a 1 board system:
% // stx (or SOP)    0x00000800
% // Packet ID       0x00000001
% // siz             0x0000002A   (# A/D data ints + #extra wds = 32+10)
% // TimeStamp Hi    1 32 bit word
% // TimeStamp Low   1 32 bit word
% // cpu             1 32 bit word
% // ttl             1 32 bit word
% // 10 extras      10 32 bit words
% // A/D data       32 32 bit words per board
% // CRC             1 32 bit XOR of the entire packet including stx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

needhdr = (nargin==1);
needdat = (nargin>=2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% since Cheetah 4.80 the DMS log files have a 16384 byte header
% determine whether that is the case in this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename, 'rb', 'ieee-le');
buf = fread(fid, 10, 'uchar=>uchar');
if all(buf(1:4)==35)
  % the file starts with "####"
  hdroffset = 16384;
elseif all(buf>31 & buf<127)
  % the file does not start with "####" but the header seems to be ascii anyway
  hdroffset = 16384;
else
  hdroffset = 0;
end

if hdroffset
  % read the ascii header from the file, it does not contain much usefull information
  orig = neuralynx_getheader(filename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% since Cheetah 5.0 the DMA log file does not always start with a complete
% packet but can start with an incomplete "junk" packet after the header.
% The junk due to that incomplete packet should then also be skipped.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hdroffset
  % skip the standard ascii header
  fseek(fid, hdroffset, 'bof');
  % read a single block
  buf = fread(fid, 274, 'uint32=>uint32');
  % search for the first known values of a proper packet
  junk = find((buf(1:(end-1))==2048) & (buf(2:(end-0))==1), 1) - 1;
  if ~isempty(junk)
    % add the additional junk samples to the header offset
    hdroffset  = hdroffset + junk*4;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% since modifying our hardware from 256 to 32 channel (19 March 2007), the
% number of channels in a DMA log file can vary. Start with the assumption
% that the file is from the 256 channel system with 8 boards that we have
% in Nijmegen and verify the CRC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nboards    = 8;
nchanboard = 32;
blocksize  = nboards*nchanboard+18; % in 4 byte words
fseek(fid, hdroffset, 'bof');
% read a single block
buf = fread(fid, blocksize, 'uint32=>uint32');
% determine the number of boards by systematically checking the CRC value
while (nboards>0)
  if neuralynx_crc(buf(1:blocksize-1))==buf(blocksize)
    % the CRC is correct for this number of boards
    break
  else
    nboards    = nboards - 1;           % decrease by one
    blocksize  = nboards*nchanboard+18; % in 4 byte words
  end
end
fclose(fid);

if ~nboards
  ft_error('could not determine the number of boards and channels');
end

% deal with ascii and numeric input for the channel selection
if nargin>3 && ~ischar(channel)
  chanindx = channel;
  channel = 'selection';
elseif nargin>3 && ischar(channel)
  switch channel
    case 'stx'
      chanindx = 1;
    case 'pid'
      chanindx = 2;
    case 'siz'
      chanindx = 3;
    case 'tsh'
      chanindx = 4;
    case 'tsl'
      chanindx = 5;
    case 'cpu'
      chanindx = 6;
    case 'ttl'
      chanindx = 7;
    case 'x01'
      chanindx = 8;
    case 'x02'
      chanindx = 9;
    case 'x03'
      chanindx = 10;
    case 'x04'
      chanindx = 11;
    case 'x05'
      chanindx = 12;
    case 'x06'
      chanindx = 13;
    case 'x07'
      chanindx = 14;
    case 'x08'
      chanindx = 15;
    case 'x09'
      chanindx = 16;
    case 'x10'
      chanindx = 17;
    case 'crc'
      chanindx = blocksize;
    case 'all'
      chanindx = 1:blocksize;
    otherwise
      ft_error('unknown value in channel');
  end
elseif nargin<4
  channel = 'all';
  chanindx = 1:blocksize;
end

if needhdr
  % these have to be hardcoded, since the DMA logfile does not contain header information
  hdr              = [];
  if hdroffset
    % remember the ascii header from the file, it does not contain much usefull information
    hdr.orig       = orig;
  end
  hdr.Fs           = 32556;                   % sampling frequency
  hdr.nChans       = nboards*nchanboard;      % number of analog channels
  hdr.nSamples     = inf;                     % number of samples per trial
  hdr.nSamplesPre  = 0;                       % number of pre-trigger samples in each trial
  hdr.nTrials      = 1;                       % number of trials
  hdr.label        = {};                      % cell-array with labels of each channel
  for i=1:hdr.nChans
    chanlabel{i} = sprintf('csc%03d', i);
  end
  statuslabel = {
    'stx'
    'pid'
    'siz'
    'tsh'
    'tsl'
    'cpu'
    'ttl'
    'x01'
    'x02'
    'x03'
    'x04'
    'x05'
    'x06'
    'x07'
    'x08'
    'x09'
    'x10'
    };
  hdr.label       = cat(1, statuslabel(:), chanlabel(:), {'crc'});  % concatenate all channel labels
  hdr.nChans      = length(hdr.label);                              % this includes all channels
  % determine the length of the file, expressed in samples
  fid = fopen(filename, 'rb', 'ieee-le');
  fseek(fid, 0, 'eof');
  hdr.nSamples = floor((ftell(fid)-hdroffset)/(blocksize*4));

  % determine the timestamp of the first and of the last sample
  datoffset = 0*blocksize;
  fseek(fid, hdroffset + datoffset, 'bof');
  buf = fread(fid, blocksize, 'uint32=>uint32');
  beg_stx = buf(1);
  beg_tsh = buf(4);
  beg_tsl = buf(5);
  datoffset = (hdr.nSamples-1)*blocksize*4;
  fseek(fid, hdroffset + datoffset, 'bof');
  buf = fread(fid, blocksize, 'uint32=>uint32');
  end_stx = buf(1);
  end_tsh = buf(4);
  end_tsl = buf(5);
  fclose(fid);
  % check that the junk at the beginning was correctly detected
  if (beg_stx~=2048)
    ft_error('problem with STX at the begin of the file');
  end
  % check that the file is truely continuous, i.e. no gaps with junk in between
  if (end_stx~=2048)
    ft_error('problem with STX at the end of the file');
  end
  % combine the two uint32 words into a uint64
  hdr.FirstTimeStamp = timestamp_neuralynx(beg_tsl, beg_tsh);
  hdr.LastTimeStamp  = timestamp_neuralynx(end_tsl, end_tsh);
  hdr.TimeStampPerSample = double(hdr.LastTimeStamp-hdr.FirstTimeStamp)./(hdr.nSamples-1);  % this should be double, since it can be fractional

  % only return the header information
  dat = hdr;
  dat.orig.Offset = hdroffset;

elseif length(chanindx)>1
  % read the desired samples for all channels
  fid = fopen(filename, 'rb', 'ieee-le');
  datoffset = (begsample-1)*blocksize*4;
  fseek(fid, hdroffset + datoffset, 'bof');  % skip the header and samples that do not have to be read
  % read the data
  dat = fread(fid, [blocksize (endsample-begsample+1)], 'int32=>int32');
  fclose(fid);
  if size(dat,2)<(endsample-begsample+1)
    ft_error('could not read all samples');
  end
  if ~strcmp(channel, 'all')
    % select the subset of desired channels
    dat = dat(chanindx,:);
  end

elseif length(chanindx)==1
  % read the desired samples for only one channel
  fid = fopen(filename, 'rb', 'ieee-le');
  datoffset = (begsample-1)*blocksize*4;
  fseek(fid, hdroffset + datoffset, 'bof');  % skip the header and samples that do not have to be read
  fseek(fid, (chanindx-1)*4, 'cof');         % skip to the channel of interest
  % Note that the last block with 274*4 bytes can sometimes be incomplete, which is relevant when endsample=inf
  dat = fread(fid, [1 (endsample-begsample+1)], 'int32=>int32', (nboards*32+18-1)*4);
  if size(dat,2)<(endsample-begsample+1)
    ft_error('could not read all samples');
  end
  fclose(fid);
end

