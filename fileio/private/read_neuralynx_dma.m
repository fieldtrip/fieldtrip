function [dat] = read_neuralynx_dma(filename, begsample, endsample, channel);

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
% $Log: read_neuralynx_dma.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.21  2008/12/15 09:40:33  roboos
% added comment in help about scaling in AD values
%
% Revision 1.20  2008/06/13 12:47:23  roboos
% added check on STX to detect discontinuous recordings (and give error)
%
% Revision 1.19  2007/12/12 11:30:28  roboos
% remember the offset of the data (header+jumk) in hdr.orig.Offset
%
% Revision 1.18  2007/10/31 14:11:26  roboos
% fixed CRC problem with the file format that was introduced in Cheetah 5.0, where the first packet could be incomplete and hence all subsequent packets are shifted
% fixed bug in getting the LastTimestamp (read as uint32 instead of int32)
%
% Revision 1.17  2007/10/25 12:48:27  roboos
% moved the reading of the ascii header to another place in the code, no functional change
%
% Revision 1.16  2007/06/13 18:53:07  roboos
% fixed bug: in the detection of header size and channel count the file was opened but not closed due to missing fclose()
%
% Revision 1.15  2007/04/02 15:12:56  roboos
% fixed detextion of ascii header (>32 should be >31)
%
% Revision 1.14  2007/03/26 12:38:19  roboos
% also detect ascii headers that do not start with '####'
%
% Revision 1.13  2007/03/19 16:57:36  roboos
% fixed bug due to hard-coded number of boards (nchans is now flexible and not always 274)
%
% Revision 1.12  2007/02/21 09:55:06  roboos
% fixed bug in TimeStampsPerSample (should be divided by nsamples-1)
%
% Revision 1.11  2007/02/20 08:55:00  roboos
% implemented automatic detection of the number of boards and channels by using the CRC value
% added the ascii header to the output
%
% Revision 1.10  2007/01/09 09:41:49  roboos
% read low and high timestamp as uint32, use seperate function to combaine them into a uint64
%
% Revision 1.9  2006/12/12 11:37:08  roboos
% read either header or data from the file
% treat datachannels and special channels the same (i.e. a file now consists of 274 'channels')
% cleaned up code, made more consistent with other functions
% allow specification of channels using either channel labels, or using channel indices
%
% Revision 1.8  2006/04/24 12:19:19  roboos
% added support for reading all channels, including the status channels
%
% Revision 1.7  2006/03/14 10:28:45  roboos
% fixed bug in read_neuralynx_dma: when reading a single (special) channel with begsample~=1 it did not skip to the begin sample
%
% Revision 1.6  2006/03/10 12:26:08  roboos
% fixed bug in reading of LastTimeStamp
%
% Revision 1.5  2006/03/06 09:50:45  roboos
% added support for reading a single special channel at a time, using the skip parameter in fread (usefull for ttl)
%
% Revision 1.4  2006/03/02 08:46:11  roboos
% account for the optional 16384 byte header introduced in v4.80
%
% Revision 1.3  2005/10/14 07:05:10  roboos
% added the first and last timestamp to the output header
%
% Revision 1.2  2005/09/09 12:32:00  roboos
% changed spacing and indentation, no functional change
%
% Revision 1.1  2005/09/05 13:05:44  roboos
% new implementation

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
  error('could not determine the number of boards and channels');
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
      error('unknown value in channel');
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
    error('problem with STX at the begin of the file');
  end
  % check that the file is truely continuous, i.e. no gaps with junk in between
  if (end_stx~=2048)
    error('problem with STX at the end of the file');
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
    error('could not read all samples');
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
    error('could not read all samples');
  end
  fclose(fid);
end

