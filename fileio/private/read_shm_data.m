function [dat, dimord] = read_shm_data(hdr, chanindx, begtrial, endtrial)

% READ_SHM_DATA reads the data in real-time from shared memory
% this is a helper function for READ_DATA

% Copyright (C) 2007-2009, Robert Oostenveld
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

% this persistent variable is used for caching the 600 packets
% which inproves the throughput when reading overlapping data segments
persistent ctf_shm

% read the data from shared memory, first the meta information only
[msgType msgId sampleNumber numSamples numChannels] = read_ctf_shm;

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

% allocate memory for the data, fill with NaNs
dat = nan(length(chanindx), hdr.nSamples, endtrial-begtrial+1);

% determine which trials/packets to read   
sel = find((trlNum>=begtrial) & (trlNum<=endtrial));

% this is for calibrating the integer values to physical values
if isfield(hdr, 'orig') && isfield(hdr.orig, 'gainV')
  gain = diag(hdr.orig.gainV(chanindx));
elseif isfield(hdr, 'orig') && isfield(hdr.orig, 'res4')
  gain = diag([hdr.orig.res4.senres(chanindx).gain]);
end

if length(chanindx)>1
  % this speeds up the multiplication, but would result in a sparse matrix when nchans=1
  gain = sparse(gain);
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
