function [hdr] = read_shm_header(filename)

% READ_SHM_HEADER reads the header in real-time from shared memory
% this is a helper function for READ_HEADER

% Copyright (C) 2007-2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% these are for remembering the header details to speed up subsequent calls
persistent previous_headerfile previous_header

% this global variable is used for caching in read_data, to improve the throughput when reading overlapping data segments
global ctf_shm
ctf_shm = [];

% decode the filename, which looks like shm://<filename>
headerfile = filetype_check_uri(filename);

if ~isempty(headerfile)
  % the headerfile has been specified by the user
  buf = [];
else
  % get the name of the headerfile from shared memory
  [msgType msgId sampleNumber numSamples numChannels] = read_ctf_shm;
  sel = find(msgType==0);
  if isempty(sel)
    error('could not determine header file location from shared memory');
  end
  buf = read_ctf_shm(sel);
  str = char(typecast(buf, 'uint8'));
  pad = find(str==0);
  headerfile = char(str(1:(pad(1)-1)));
end

if isempty(previous_headerfile) || isempty(previous_header)

  % read the header details and remember for the subsequent calls
  hdr = read_header(headerfile, 'cache', true);
  previous_headerfile = headerfile;
  previous_header     = hdr;

  % inform AcqBuffer about the trigger channels, to allow flank detection
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

else
  % don't read the header details again, only update the number of samples/trials
  hdr = previous_header;
end

% the following information is determined from shared memory
sel = find(msgType==1);  % these are the data packets
hdr.nTrials     = double(max(sampleNumber(sel)))/double(max(numSamples(sel)));
hdr.nSamples    = double(max(numSamples(sel)));
hdr.nSamplesPre = 0;

