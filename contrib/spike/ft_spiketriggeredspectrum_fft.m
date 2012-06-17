function [sts] = ft_spiketriggeredspectrum_fft(cfg, data)

% FT_SPIKETRIGGEREDSPECTRUM_FFT computes the Fourier spectrup of the LFP around
% the spikes.
%
% Use as
%   [sts] = ft_spiketriggeredspectrum_fft(cfg, data)
%
% The input data should be organised in a structure as obtained from
% the FT_APPENDSPIKE. The configuration should be according to
%
%   cfg.timwin       = [begin end], time around each spike (default = [-0.1 0.1])
%   cfg.foilim       = [begin end], frequency band of interest (default = [0 150])
%   cfg.taper        = 'dpss', 'hanning' or many others, see WINDOW (default = 'hanning')
%   cfg.tapsmofrq    = number, the amount of spectral smoothing through
%                      multi-tapering. Note that 4 Hz smoothing means
%                      plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
%                      Note: multitapering rotates phases (no problem for consistency)
%   cfg.spikechannel = string, name of single spike channel to trigger on
%   cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                      see FT_CHANNELSELECTION for details
%   cfg.feedback     = 'no', 'text', 'textbar', 'gui' (default = 'no')
%
% If the triggered spike leads a spike in another channel, then the angle
% of the Fourier spectrum of that other channel will be negative. 
% Earlier phases are in clockwise direction.
% See FT_SPIKETRIGGEREDINTERPOLATION to remove segments of LFP around
% spikes
% See FT_SPIKETRIGGEREDSPECTRUM_CONVOL for an alternative implementation
% based on convolution
%
% Output sts can be input to FT_SPIKETRIGGEREDSPECTRUM_STAT
%
% A phase of zero corresponds to the spike being on the peak of the LFP
% oscillation.
% A phase of 180 degree corresponds to the spike being in the through of the
% oscillation.
% A phase of 45 degrees corresponds to the spike being just after the
% peak in the LFP.
% This function uses a NaN-aware spectral estimation technique, which will
% default to the standard Matlab FFT routine if no NaNs are present. The
% fft_along_rows subfunction below demonstrates the expected function
% behaviour.

% Copyright (C) 2008, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% check input data structure
data = ft_checkdata(data,'datatype', 'raw', 'feedback', 'yes');

% these were supported in the past, but are not any more (for consistency with other spike functions)
cfg = ft_checkconfig(cfg, 'forbidden', {'inputfile','outputfile'}); 

%get the options
cfg.timwin       = ft_getopt(cfg, 'timwin',[-0.1 0.1]);
cfg.spikechannel = ft_getopt(cfg,'spikechannel', []);
cfg.channel      = ft_getopt(cfg,'channel', 'all');
cfg.feedback     = ft_getopt(cfg,'feedback', 'no');
cfg.tapsmofrq    = ft_getopt(cfg,'tapsmofrq', 4);
cfg.taper        = ft_getopt(cfg,'taper', 'hanning');
cfg.foilim       = ft_getopt(cfg,'foilim', [0 150]);

% ensure that the options are valid
cfg = ft_checkopt(cfg,'timwin','doublevector');
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double', 'empty'});
cfg = ft_checkopt(cfg,'channel', {'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'feedback', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'taper', 'char');
cfg = ft_checkopt(cfg,'tapsmofrq', 'doublescalar');
cfg = ft_checkopt(cfg,'foilim', 'doublevector');

if strcmp(cfg.taper, 'sine')
  error('sorry, sine taper is not yet implemented');
end

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
cfg.channel = ft_channelselection(cfg.channel, data.label);
chansel     = match_str(data.label, cfg.channel);
nchansel    = length(cfg.channel);  % number of channels

% determine the spike channel on which will be triggered
cfg.spikechannel = ft_channelselection(cfg.spikechannel, data.label);
spikesel         = match_str(data.label, cfg.spikechannel);
nspikesel        = length(cfg.spikechannel);    % number of channels

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
if ~strcmp(cfg.taper,'dpss')
  taper  = window(cfg.taper, numsmp);
  taper  = taper./norm(taper);
else
  % not implemented yet: keep tapers, or selecting only a subset of them.
  taper  = dpss(numsmp, cfg.tapsmofrq);
  taper  = taper(:,1:end-1);            % we get 2*NW-1 tapers
  taper  = sum(taper,2)./size(taper,2); % using the linearity of multitapering
end
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
time  = randn(1,numsmp); % this is actually not used
spike(1-begpad) = 1;
spike_fft = specest_nanfft(spike, time);
spike_fft = spike_fft(fbeg:fend);
spike_fft = spike_fft./abs(spike_fft);
rephase   = sparse(diag(conj(spike_fft)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ft_progress('init', 'text',     'Please wait...');
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
  
  spectrum{i} = zeros(length(spikesmp), nchansel, fend-fbeg+1);
  ft_progress(i/ntrial, 'spectrally decomposing data for trial %d of %d, %d spikes', i, ntrial, length(spikesmp));  
  for j=1:length(spikesmp)
    begsmp = spikesmp(j) + begpad;
    endsmp = spikesmp(j) + endpad;
    
    if (begsmp<1)
      segment = nan*zeros(nchansel, numsmp);
    elseif endsmp>size(data.trial{i},2)
      segment = nan*zeros(nchansel, numsmp);
    else
      segment = data.trial{i}(chansel,begsmp:endsmp);
    end
    
    % substract the DC component from every segment, to avoid any leakage of the taper
    segmentMean = repmat(nanmean(segment,2),1,numsmp); % nChan x Numsmp
    segment     = segment - segmentMean; % LFP has average of zero now (no DC)
    
    time  = randn(size(segment)); % this is actually not used
    
    % taper the data segment around the spike and compute the fft
    segment_fft = specest_nanfft(segment * taper, time);
    
    % select the desired output frquencies and normalize
    segment_fft = segment_fft(:,fbeg:fend) ./ sqrt(numsmp/2);
    
    % rotate the estimated phase at each frequency to correct for the segment t=0 not being at the first sample
    segment_fft = segment_fft * rephase;
    
    % store the result for this spike in this trial
    spectrum{i}(j,:,:) = segment_fft;
    
  end % for each spike in this trial  
end % for each trial
ft_progress('close');
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sts.lfplabel       = data.label(chansel);
sts.freq           = freqaxis(fbeg:fend);
sts.dimord         = 'rpt_chan_freq';
sts.fourierspctrm  = {cat(1, spectrum{:})};
sts.time           = {cat(2,spiketime{:})'};  % this deviates from the standard output, but is included for reference
sts.trial          = {cat(2,spiketrial{:})'}; % this deviates from the standard output, but is included for reference
sts.dimord = '{chan}_spike_lfpchan_freq';
for i = 1:ntrial
  sts.trialtime(i,:) = [data.time{i}(1) data.time{i}(end)];
end
sts.label   = data.label(spikesel);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous data
ft_postamble history sts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for demonstration purpose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = fft_along_rows(x)
y = fft(x, [], 2); % use normal Matlab function to compute the fft along 2nd dimension

