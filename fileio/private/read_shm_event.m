function [event] = read_shm_event(filename, varargin);

% READ_SHM_EVENT reads the events in real-time from shared memory
% this is a helper function for READ_EVENT

% Copyright (C) 2007, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% get the optional input arguments
hdr       = keyval('header',    varargin);
type      = keyval('type',      varargin);
minsample = keyval('minsample', varargin);
maxsample = keyval('maxsample', varargin);

if isempty(hdr)
  hdr = read_header(filename);
end

% Acq is writing the data to shared memory in real-time
% here we read the data from shared memory, first the meta information only
[msgType msgId sampleNumber numSamples numChannels] = read_ctf_shm;

event = [];

% generate 'trial' events for each packet
if isempty(type) || any(strcmp('trial', type))
  % there seems to be a bug in Acq, causing the messageId to wrap around
  % hence it cannot be used as index into the packets, so construct a new trial numbering vector
  trlNum = nan(size(msgId));
  trlNum(msgType==1) = sampleNumber(msgType==1)./numSamples(msgType==1);
  % make an event for each data packet
  sel = find(msgType==1);
  for i=1:length(sel)
    event(end+1).type = 'trial';
    event(end ).value    = [];
    event(end ).sample   = double(sampleNumber(sel(i))+1);  % offset by one
    event(end ).duration = double(numSamples(sel(i)));
    event(end ).offset   = 0;
  end
end 

if any(msgType==0)
  % there is a setup packet, hence assume that realtime trigger detection has been enabled
  % read the triggers that were detected by AcqBuffer
  sel = find(msgType==0);
  if ~isempty(sel)
    buf = read_ctf_shm(sel);
    % the last 9000 samples in the setup buffer contain up to 3x3000 trigger values
    % the first row contains the channel, the second row the sample and the tird row the value
    trg = reshape(buf((end-9000+1):end), 3, 3000);
    trg = trg(:, trg(2,:)>0);
    % apply some event filters
    if ~isempty(minsample)
      trg = trg(:, trg(2,:)>=minsample);
    end
    if ~isempty(maxsample)
      trg = trg(:, trg(2,:)<=maxsample);
    end
    for i=1:size(trg,2)
      event(end+1).type   = hdr.label{trg(1,i)+1};  % zero-offset
      event(end  ).sample = trg(2,i)+1;             % zero-offset
      event(end  ).value  = trg(3,i);
    end
  end
else
  % determine the trigger channels from the header
  if isfield(hdr, 'orig') && isfield(hdr.orig, 'sensType')
    % determine the samples that are available in shared memory
    sel = find(msgType==1);
    begsample = min(sampleNumber(sel))+1;
    endsample = max(sampleNumber(sel))+hdr.nSamples;
    trgchan = find(hdr.orig.sensType(:)'==11);
    if ~isempty(type)
      % only look in the selected trigger channels
      trgchan = match_str(hdr.label, intersect(hdr.label(trgchan), type));
    end
    for i=trgchan(:)'
      % read this trigger channel as raw data, can safely assume that it is continuous
      trig = read_data(filename, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', i, 'checkboundary', false);
      % correct for reading it as signed integer, whereas it should be an unsigned int
      trig(find(trig<0)) = trig(find(trig<0)) + 2^32;
      % TODO apply some event filters
      % convert the trigger into an event with a value at a specific sample
      for j=find(diff([0 trig(:)'])>0)
        event(end+1).type   = hdr.label{i};
        event(end  ).value  = trig(j);
        event(end  ).sample = j+begsample-1;
      end
    end
  end
end

