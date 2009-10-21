function [freq] = spiketriggeredspectrum(cfg, data)

% SPIKETRIGGEREDSPECTRUM computes the Fourier spectrup of the LFP around the spikes.
%
% Use as
%   [freq] = spiketriggeredspectrum(cfg, data)
%
% The input data should be organised in a structure as obtained from
% the APPENDSPIKE function. The configuration should be according to
%
%   cfg.timwin       = [begin end], time around each spike (default = [-0.1 0.1])
%   cfg.foilim       = [begin end], frequency band of interest (default = [0 150])
%   cfg.taper        = 'dpss', 'hanning' or many others, see WINDOW (default = 'hanning')
%   cfg.spikechannel = string, name of single spike channel to trigger on
%   cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                      see CHANNELSELECTION for details
%   cfg.feedback     = 'no', 'text', 'textbar', 'gui' (default = 'no')
%
% If the triggered spike leads a spike in another channel, then the angle
% of the Fourier spectrum of that other channel will be negative. NOTE that
% this should be checked for consistency.

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: spiketriggeredspectrum.m,v $
% Revision 1.8  2009/01/20 16:53:29  roboos
% fixed problem with fft handle not expecting key-value input arguments (thanks to Thilo)
%
% Revision 1.7  2008/10/01 08:23:05  roboos
% use specest_nanfft as default
%
% Revision 1.6  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.5  2008/07/15 18:17:57  roboos
% fixed documentation
%
% Revision 1.4  2008/04/09 14:39:22  roboos
% remove trials without data (too close to trial border) from the output
%
% Revision 1.3  2008/04/09 14:23:15  roboos
% replicate double spikes, for keeptrials
% give error if spike count >5
%
% Revision 1.2  2008/03/18 21:58:18  roboos
% added normalization of the spectrum by the number of samples, and rotation of the phase to correct for t=0 being in the middle of the segment instead of at the begin
%
% Revision 1.1  2008/03/17 16:42:35  roboos
% initial implementation, phase rotation still needs to be built in
%

fieldtripdefs

% set the defaults
if ~isfield(cfg, 'timwin'),       cfg.timwin = [-0.1 0.1];    end
if ~isfield(cfg, 'foilim'),       cfg.foilim = [0 150];       end
if ~isfield(cfg, 'taper'),        cfg.taper = 'hanning';      end
if ~isfield(cfg, 'channel'),      cfg.channel = 'all';        end
if ~isfield(cfg, 'spikechannel'), cfg.spikechannel = [];      end
% if ~isfield(cfg, 'keeptrials'),   cfg.keeptrials = 'no';      end
if ~isfield(cfg, 'feedback'),     cfg.feedback = 'no';        end

if strcmp(cfg.taper, 'dpss') && strcmp(cfg.taper, 'sine')
  error('sorry, multitapering is not yet implemented');
end

% see subfunction below, which demonstrates the expected function behaviour
% my_fft = @fft_along_rows; 

% use a NaN-aware spectral estimation technique, which will default to the
% standard Matlab FFT routine if no NaNs are present
my_fft = @specest_nanfft;

% autodetect the spike channels
ntrial = length(data.trial);
nchans  = length(data.label);
spikechan = zeros(nchans,1);
for i=1:ntrial
  for j=1:nchans
    spikechan(j) = spikechan(j) + all(data.trial{i}(j,:)==0 | data.trial{i}(j,:)==1 | data.trial{i}(j,:)==2);
  end
end
spikechan = (spikechan==ntrial);

% determine the channels to be averaged
cfg.channel = channelselection(cfg.channel, data.label);
chansel     = match_str(data.label, cfg.channel);
nchansel    = length(cfg.channel);	% number of channels

% determine the spike channel on which will be triggered
cfg.spikechannel = channelselection(cfg.spikechannel, data.label);
spikesel         = match_str(data.label, cfg.spikechannel);
nspikesel        = length(cfg.spikechannel);	% number of channels

if nspikesel==0
  error('no spike channel selected');
end

if nspikesel>1
  error('only supported for a single spike channel');
end

if ~spikechan(spikesel)
  error('the selected spike channel seems to contain continuous data');
end

begpad = round(cfg.timwin(1)*data.fsample);
endpad = round(cfg.timwin(2)*data.fsample);
numsmp = endpad - begpad + 1;
taper  = window('hanning', numsmp);
taper  = taper./norm(taper);
taper  = sparse(diag(taper));

spectrum    = cell(1,ntrial);
spiketime   = cell(1,ntrial);
spiketrial  = cell(1,ntrial);
cumsum = zeros(nchansel, numsmp);
cumcnt = 0;

timeaxis = linspace(cfg.timwin(1),cfg.timwin(2), numsmp);
freqaxis = linspace(0, data.fsample, numsmp);
fbeg = nearest(freqaxis, cfg.foilim(1));
fend = nearest(freqaxis, cfg.foilim(2));
% update the configuration to acocunt for rounding off differences
cfg.foilim(1) = freqaxis(fbeg);
cfg.foilim(2) = freqaxis(fend);

% make a representation of the spike, this is used for the phase rotation
spike = zeros(1,numsmp);
spike(1-begpad) = 1;
spike_fft = my_fft(spike);
spike_fft = spike_fft(fbeg:fend);
spike_fft = spike_fft./abs(spike_fft);
rephase   = sparse(diag(conj(spike_fft)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:ntrial
  spikesmp = find(data.trial{i}(spikesel,:));
  spikecnt = data.trial{i}(spikesel,spikesmp);

  if any(spikecnt>5) || any(spikecnt<0)
    error('the spike count lies out of the regular bounds');
  end

  % instead of doing the bookkeeping of double spikes below, replicate the double spikes by looking at spikecnt
  sel = find(spikecnt>1);
  tmp = zeros(1,sum(spikecnt(sel)));
  n   = 1;
  for j=1:length(sel)
    for k=1:spikecnt(sel(j))
      tmp(n) = spikesmp(sel(j));
      n = n + 1;
    end
  end
  spikesmp(sel) = [];                     % remove the double spikes
  spikecnt(sel) = [];                     % remove the double spikes
  spikesmp = [spikesmp tmp];              % add the double spikes as replicated single spikes
  spikecnt = [spikecnt ones(size(tmp))];  % add the double spikes as replicated single spikes
  spikesmp = sort(spikesmp);              % sort them to keep the original ordering (not needed on spikecnt, since that is all ones)

  spiketime{i}  = data.time{i}(spikesmp);
  spiketrial{i} = i*ones(size(spikesmp));
  fprintf('processing trial %d of %d (%d spikes)\n', i, ntrial, sum(spikecnt));

  spectrum{i} = zeros(length(spikesmp), nchansel, fend-fbeg+1);

  progress('init', cfg.feedback, 'spectrally decomposing data around spikes');
  for j=1:length(spikesmp)
    progress(i/ntrial, 'spectrally decomposing data around spike %d of %d\n', j, length(spikesmp));
    begsmp = spikesmp(j) + begpad;
    endsmp = spikesmp(j) + endpad;

    if (begsmp<1)
      segment = nan*zeros(nchansel, numsmp);
    elseif endsmp>size(data.trial{i},2)
      segment = nan*zeros(nchansel, numsmp);
    else
      segment = data.trial{i}(chansel,begsmp:endsmp);
    end

    % taper the data segment around the spike and compute the fft
    segment_fft = my_fft(segment * taper);

    % select the desired output frquencies and normalize
    segment_fft = segment_fft(:,fbeg:fend) ./ sqrt(numsmp/2);

    % rotate the estimated phase at each frequency to correct for the segment t=0 not being at the first sample
    segment_fft = segment_fft * rephase;

    % store the result for this spike in this trial
    spectrum{i}(j,:,:) = segment_fft;

  end % for each spike in this trial
  progress('close');

end % for each trial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq.label          = data.label(chansel);
freq.freq           = freqaxis(fbeg:fend);
freq.dimord         = 'rpt_chan_freq';

freq.fourierspctrm  = cat(1, spectrum{:});
freq.origtime       = cat(2,spiketime{:})';  % this deviates from the standard output, but is included for reference
freq.origtrial      = cat(2,spiketrial{:})'; % this deviates from the standard output, but is included for reference

% select all trials that do not contain data in the first sample
sel = isnan(freq.fourierspctrm(:,1,1));
fprintf('removing %d trials from the output that do not contain data\n', sum(sel));
% remove the selected trials from the output
freq.fourierspctrm  = freq.fourierspctrm(~sel,:,:);
freq.origtime       = freq.origtime(~sel);
freq.origtrial      = freq.origtrial(~sel);

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: spiketriggeredspectrum.m,v 1.8 2009/01/20 16:53:29 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
freq.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = fft_along_rows(x)
y = fft(x, [], 2); % use normal Matlab function to compute the fft along 2nd dimension 

