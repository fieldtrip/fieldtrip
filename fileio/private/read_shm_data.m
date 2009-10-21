function [dat, dimord] = read_shm_data(hdr, chanindx, begtrial, endtrial)

% READ_SHM_DATA reads the data in real-time from shared memory
% this is a helper function for READ_DATA

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: read_shm_data.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.2  2008/11/20 13:56:11  roboos
% fixed gain when header is read using ctf p-files
%
% Revision 1.1  2007/08/01 12:12:06  roboos
% moved the actual code from the normal functions into these helper functions
% fixed some bugs related to reading the header from a user-specified res4 file and setting the trigger detection
% use the new trigger detection in AcqBuffer when possible, use the old trigger detection if no setup buffer is present
% implemented caching of the data packets (using global variable ctf_shm)
%

% read the data from shared memory, first the meta information only
[msgType msgId sampleNumber numSamples numChannels] = read_ctf_shm;

% this global variable is used for caching in read_data
% which inproved the throughput when reading overlapping data segments
global ctf_shm
if isempty(ctf_shm)
  ctf_shm.msgType      = nan(size(msgType));
  ctf_shm.msgId        = nan(size(msgId));
  ctf_shm.sampleNumber = nan(size(sampleNumber));
  ctf_shm.numSamples   = nan(size(numSamples));
  ctf_shm.numChannels  = nan(size(numChannels));
  ctf_shm.data = {};
  ctf_shm.hit   = 0;
  ctf_shm.mis   = 0;
end

% there seems to be a bug in Acq, causing the messageId to wrap around
% hence it cannot be used as index into the packets, so construct a new trial numbering vector 
trlNum = nan(size(msgId));
trlNum(msgType==1) = sampleNumber(msgType==1)./numSamples(msgType==1) + 1;

% allocate memory for teh data, fill with NaNs
dat = nan(length(chanindx), hdr.nSamples, endtrial-begtrial+1);

% determine which trials/packets to read   
sel = find((trlNum>=begtrial) & (trlNum<=endtrial));

% this is for calibrating the integer values to physical values
if isfield(hdr, 'orig') && isfield(hdr.orig, 'gainV')
  gain = sparse(diag(hdr.orig.gainV(chanindx)));
elseif isfield(hdr, 'orig') && isfield(hdr.orig, 'res4')
  gain = sparse(diag([hdr.orig.res4.senres(chanindx).gain]));
end
  
for i=1:length(sel)
  % nchan = numChannels(i);
  % nsmp  = numSamples(i);
  nchan = hdr.nChans;
  nsmp  = hdr.nSamples;
  % buf = read_ctf_shm(sel(i));
  % buf = buf(1:nchan*nsmp);
  if all([msgType(sel(i)) msgId(sel(i)) sampleNumber(sel(i)) numSamples(sel(i)) numChannels(sel(i))] ==  [ctf_shm.msgType(sel(i)) ctf_shm.msgId(sel(i)) ctf_shm.sampleNumber(sel(i)) ctf_shm.numSamples(sel(i)) ctf_shm.numChannels(sel(i))])
    % get the data from the cache
    buf = ctf_shm.data{sel(i)};
    ctf_shm.hit = ctf_shm.hit+1;
  else
    % read the data from shared memory, update the cache
    buf = read_ctf_shm(sel(i),nchan*nsmp);
    ctf_shm.data{sel(i)}         = buf;
    ctf_shm.msgType(sel(i))      = msgType(sel(i));
    ctf_shm.msgId(sel(i))        = msgId(sel(i));
    ctf_shm.sampleNumber(sel(i)) = sampleNumber(sel(i));
    ctf_shm.numSamples(sel(i))   = numSamples(sel(i));
    ctf_shm.numChannels(sel(i))  = numChannels(sel(i));
    ctf_shm.mis                  = ctf_shm.mis+1;
  end
  buf = reshape(buf, nchan, nsmp);
  thistrial = trlNum(sel(i));
  % apply calibration to the selected channels
  dat(:,:,thistrial-begtrial+1) = gain * double(buf(chanindx,:));
end

% if any(isnan(dat(:)))
%   warning('data has been padded with NaNs');
%   fprintf('trials present   = %d - %d\n', min(trlNum), max(trlNum));
%   fprintf('trials requested = %d - %d\n', begtrial, endtrial);
% end

dimord = 'chans_samples_trials';
