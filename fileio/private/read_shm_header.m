function [hdr] = read_shm_header(filename)

% READ_SHM_HEADER reads the header information in real-time from shared memory
% this is a helper function for READ_HEADER

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: read_shm_header.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.5  2009/01/12 12:01:55  roboos
% determine number of samples and trials from only the data blocks, not all blocks. It seems that memcpy has changed, which now sometimes seems to cause a block to be "corupt" for a small amount of time.
%
% Revision 1.4  2008/12/02 08:28:40  roboos
% fixed typo, missing )
%
% Revision 1.3  2008/11/20 13:58:11  roboos
% fixed hdr.nTrials and nSamples (ensure that it is read as blocks)
% use caching when reading header from file
%
% Revision 1.2  2008/10/08 10:08:06  roboos
% detect trigger channels also when header is read using ctf_new
%
% Revision 1.1  2007/08/01 12:12:06  roboos
% moved the actual code from the normal functions into these helper functions
% fixed some bugs related to reading the header from a user-specified res4 file and setting the trigger detection
% use the new trigger detection in AcqBuffer when possible, use the old trigger detection if no setup buffer is present
% implemented caching of the data packets (using global variable ctf_shm)
%

% this global variable is used for caching in read_data, to improve the throughput when reading overlapping data segments
global ctf_shm
ctf_shm = [];

% decode the filename, which looks like shm://<filename>
headerfile = filetype_check_uri(filename);

if ~isempty(headerfile)
  % the headerfile has been specified by the user
  hdr = read_header(headerfile, 'cache', true);
  buf = [];
  [msgType msgId sampleNumber numSamples numChannels] = read_ctf_shm;
  sel = find(msgType==0);
  if ~isempty(sel)
    buf = read_ctf_shm(sel);
  else
    buf = [];
  end
else
  % get the name of the headerfile from shared memory
  [msgType msgId sampleNumber numSamples numChannels] = read_ctf_shm;
  sel = find(msgType==0);
  if ~isempty(sel)
    buf = read_ctf_shm(sel);
    str = char(typecast(buf, 'uint8'));
    pad = find(str==0);
    headerfile = char(str(1:(pad(1)-1)));
    hdr = read_header(headerfile, 'cache', true);
  else
    error('could not determine header file location from shared memory');
  end
end

if ~isempty(buf)
  % meg channels are 5, refmag 0, refgrad 1, adcs 18, trigger 11, eeg 9
  if isfield(hdr, 'orig') && isfield(hdr.orig, 'sensType')
    origSensType = hdr.orig.sensType;
  elseif isfield(hdr, 'orig') && isfield(hdr.orig, 'res4')
    origSensType =  [hdr.orig.res4.senres.sensorTypeIndex];
  end
  % do online detection of triggers inside AcqBuffer
  trgchan = find(origSensType(:)'==11);
  if length(trgchan)>10
    error('online detection of triggers in AcqBuffer only works with up to 10 trigger channels');
  end
  % specify the channel count and the trigger channels in the setup buffer
  buf(28160-9011) = hdr.nChans;        % tell the number of actual channels
  buf(28160-9010) = length(trgchan);   % tell the number of trigger channels
  for i=1:length(trgchan)
    buf(28160-9010+i) = trgchan(i)-1;  % tell the index of the trigger channel, zero offset
  end
  % the setup packet may have moved in the meantine, determine its latest location
  [msgType msgId sampleNumber numSamples numChannels] = read_ctf_shm;
  sel = find(msgType==0);
  % write the updated setup packet, this should cause AcqBuffer to do online trigger detection
  write_ctf_shm(sel, 0, 0, 0, 0, 0, buf); 
else
  warning('no setup in shared memory, could not enable trigger detection');
end

% the following information is determined from shared memory
sel = find(msgType==1);  % these are the data packets
hdr.nTrials     = double(max(sampleNumber(sel)))/double(max(numSamples(sel)));
hdr.nSamples    = double(max(numSamples(sel)));
hdr.nSamplesPre = 0;


