function [sts] = ft_spiketriggeredspectrum_fft(cfg, data, spike)

% FT_SPIKETRIGGEREDSPECTRUM_FFT computes the Fourier spectrum (amplitude and phase)
% of the LFP around the % spikes.  A phase of zero corresponds to the spike being on
% the peak of the LFP oscillation. A phase of 180 degree corresponds to the spike being
% in the through of the oscillation. A phase of 45 degrees corresponds to the spike
% being just after the peak in the LFP.
%
% If the triggered spike leads a spike in another channel, then the angle of the Fourier
% spectrum of that other channel will be negative. Earlier phases are in clockwise
% direction. 
%
% Use as
%   [sts] = ft_spiketriggeredspectrum_convol(cfg,data,spike)
% or 
%   [sts] = ft_spiketriggeredspectrum_convol(cfg,data)
% where the spike data can either be contained in the DATA input or in the SPIKE input.
%
% The input DATA should be organised as the raw datatype, obtained from FT_PREPROCESSING
% or FT_APPENDSPIKE. 
%
% The (optional) input SPIKE should be organised as the spike or the raw datatype,
% obtained from FT_SPIKE_MAKETRIALS or FT_PREPROCESSING (in that case, conversion is done
% within the function)
%
% Important is that data.time and spike.trialtime should be referenced relative to the
% same trial trigger times.
%
% The configuration should be according to
%   cfg.timwin       = [begin end], time around each spike (default = [-0.1 0.1])
%   cfg.foilim       = [begin end], frequency band of interest (default = [0 150])
%   cfg.taper        = 'dpss', 'hanning' or many others, see WINDOW (default = 'hanning')
%   cfg.tapsmofrq    = number, the amount of spectral smoothing through
%                      multi-tapering. Note that 4 Hz smoothing means plus-minus 4 Hz,
%                      i.e. a 8 Hz smoothing box. Note: multitapering rotates phases (no
%                      problem for consistency)
%   cfg.spikechannel = string, name of spike channels to trigger on cfg.channel = Nx1
%                      cell-array with selection of channels (default = 'all'),
%                      see FT_CHANNELSELECTION for details
%   cfg.feedback     = 'no', 'text', 'textbar', 'gui' (default = 'no')
%
% The output STS data structure can be input to FT_SPIKETRIGGEREDSPECTRUM_STAT
%
% This function uses a NaN-aware spectral estimation technique, which will default to the
% standard Matlab FFT routine if no NaNs are present. The fft_along_rows subfunction below
% demonstrates the expected function behaviour.
%
% See FT_SPIKETRIGGEREDINTERPOLATION to remove segments of LFP around spikes.
% See FT_SPIKETRIGGEREDSPECTRUM_CONVOL for an alternative implementation based
% on convolution

% Copyright (C) 2008, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance data spike
ft_preamble trackconfig

% check input data structure
data = ft_checkdata(data,'datatype', 'raw', 'feedback', 'yes');
if nargin==3
  spike = ft_checkdata(spike, 'datatype', {'spike'}, 'feedback', 'yes');
end

% these were supported in the past, but are not any more (for consistency with other spike functions)
cfg = ft_checkconfig(cfg, 'forbidden', {'inputfile','outputfile'}); 

%get the options
cfg.timwin       = ft_getopt(cfg, 'timwin',[-0.1 0.1]);
cfg.spikechannel = ft_getopt(cfg,'spikechannel', 'all');
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

% get the spikechannels
if nargin==2
  
  % autodetect the spikechannels and EEG channels
  [spikechannel, eegchannel] = detectspikechan(data);
  
  % make the final selection of spike channels and check
  if strcmp(cfg.spikechannel, 'all'), 
    cfg.spikechannel = spikechannel; 
  else
    cfg.spikechannel = ft_channelselection(cfg.spikechannel, data.label);  
    if ~all(ismember(cfg.spikechannel,spikechannel)), 
      error('some selected spike channels appear eeg channels'); 
    end        
  end
  
  % make the final selection of EEG channels and check
  if strcmp(cfg.channel,'all')  
    cfg.channel = eegchannel;
  else
    cfg.channel      = ft_channelselection(cfg.channel, data.label);  
    if ~all(ismember(cfg.channel,eegchannel)), 
      warning('some of the selected eeg channels appear spike channels'); 
    end    
  end    
  
  % select the data and convert to a spike structure
  data_spk = ft_selectdata(data,'channel', cfg.spikechannel);
  data     = ft_selectdata(data,'channel', cfg.channel); % leave only LFP
  spike    = ft_checkdata(data_spk,'datatype', 'spike');
  clear data_spk % remove the continuous data
else
  cfg.spikechannel = ft_channelselection(cfg.spikechannel, spike.label);  
  cfg.channel      = ft_channelselection(cfg.channel, data.label);
end 

% determine the channel indices and number of chans
chansel          = match_str(data.label, cfg.channel); % selected channels
nchansel         = length(cfg.channel);                % number of channels
spikesel         = match_str(spike.label, cfg.spikechannel);
nspikesel        = length(spikesel); % number of spike channels
if nspikesel==0, error('no spike channel selected'); end

% construct the taper
if ~isfield(data, 'fsample'), data.fsample = 1/mean(diff(data.time{1})); end
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

% preallocate the output structures for different units / trials
ntrial      = length(data.trial);
spectrum    = cell(nspikesel,ntrial);
spiketime   = cell(nspikesel,ntrial);
spiketrial  = cell(nspikesel,ntrial);

% select the frequencies
freqaxis = linspace(0, data.fsample, numsmp);
fbeg = nearest(freqaxis, cfg.foilim(1));
fend = nearest(freqaxis, cfg.foilim(2));

% update the configuration to account for rounding off differences
cfg.foilim(1) = freqaxis(fbeg);
cfg.foilim(2) = freqaxis(fend);

% make a representation of the spike, this is used for the phase rotation
spike_repr = zeros(1,numsmp);
time       = linspace(cfg.timwin(1),cfg.timwin(2), numsmp);
spike_repr(1-begpad) = 1;
spike_fft = specest_nanfft(spike_repr, time);
spike_fft = spike_fft(fbeg:fend);
spike_fft = spike_fft./abs(spike_fft);
rephase   = sparse(diag(conj(spike_fft)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ft_progress('init', 'text',     'Please wait...');
for iUnit  = 1:nspikesel
  for iTrial = 1:ntrial

    % select the spikes that fell in the trial and convert to samples
    timeBins = [data.time{iTrial}  data.time{iTrial}(end)+1/data.fsample] - (0.5/data.fsample);      
    hasTrial = spike.trial{spikesel(iUnit)} == iTrial; % find the spikes that are in the trial
    ts       = spike.time{spikesel(iUnit)}(hasTrial); % get the spike times for these spikes
    ts       = ts(ts>=timeBins(1) & ts<=timeBins(end)); % only select those spikes that fall in the trial window
    [ignore,spikesmp] = histc(ts,timeBins);      
    if ~isempty(ts)
      ts(spikesmp==0 | spikesmp==length(timeBins)) = [];
    end
      
        
    spikesmp(spikesmp==0 | spikesmp==length(timeBins)) = [];
    
    % store in the output cell arrays as column vectors
    spiketime{iUnit, iTrial}  = ts(:);
    tr = iTrial*ones(size(spikesmp));
    spiketrial{iUnit, iTrial} = tr(:);

    % preallocate the spectrum
    spectrum{iUnit, iTrial} = zeros(length(spikesmp), nchansel, fend-fbeg+1);
    
    % compute the spiketriggered spectrum
    ft_progress(iTrial/ntrial, 'spectrally decomposing data for trial %d of %d, %d spikes for unit %d', iTrial, ntrial, length(spikesmp), iUnit);  
    for j=1:length(spikesmp)
      
      % selected samples
      begsmp = spikesmp(j) + begpad;
      endsmp = spikesmp(j) + endpad;

      % handle spikes near the borders of the trials
      if (begsmp<1)
        segment = nan(nchansel, numsmp);
      elseif endsmp>size(data.trial{iTrial},2)
        segment = nan(nchansel, numsmp);
      else
        segment = data.trial{iTrial}(chansel,begsmp:endsmp);
      end

      % substract the DC component from every segment, to avoid any leakage of the taper
      segmentMean = repmat(nanmean(segment,2),1,numsmp); % nChan x Numsmp
      segment     = segment - segmentMean; % LFP has average of zero now (no DC)

      % taper the data segment around the spike and compute the fft
      segment_fft = specest_nanfft(segment * taper, time);

      % select the desired output frquencies and normalize
      segment_fft = segment_fft(:,fbeg:fend) ./ sqrt(numsmp/2);

      % rotate the estimated phase at each frequency to correct for the segment t=0 not being at the first sample
      segment_fft = segment_fft * rephase;

      % store the result for this spike in this trial
      spectrum{iUnit, iTrial}(j,:,:) = segment_fft;

    end % for each spike in this trial  
  end % for each trial
end
ft_progress('close');
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results in a structure that is a spike structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sts.lfplabel       = data.label(chansel);
sts.freq           = freqaxis(fbeg:fend);
sts.dimord         = 'rpt_chan_freq';
for iUnit = 1:nspikesel
  sts.fourierspctrm{iUnit}  = cat(1, spectrum{iUnit,:});
  spectrum(iUnit,:) = {[]}; % free from the memory
  sts.time{iUnit}           = cat(1,spiketime{iUnit,:});  
  sts.trial{iUnit}          = cat(1,spiketrial{iUnit,:}); 
end
sts.dimord    = '{chan}_spike_lfpchan_freq';
sts.trialtime = spike.trialtime;
sts.label     = spike.label(spikesel);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance sts
ft_postamble history    sts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spikelabel, eeglabel] = detectspikechan(data)

maxRate = 1000; % default on what we still consider a neuronal signal

% autodetect the spike channels
ntrial = length(data.trial);
nchans  = length(data.label);
spikechan = zeros(nchans,1);
for i=1:ntrial
  for j=1:nchans
    hasAllInts    = all(isnan(data.trial{i}(j,:)) | data.trial{i}(j,:) == round(data.trial{i}(j,:)));
    hasAllPosInts = all(isnan(data.trial{i}(j,:)) | data.trial{i}(j,:)>=0);
    fr            = nansum(data.trial{i}(j,:),2) ./ (data.time{i}(end)-data.time{i}(1));    
    spikechan(j) = spikechan(j) + double(hasAllInts & hasAllPosInts & fr<=maxRate);
  end
end
spikechan = (spikechan==ntrial);

spikelabel = data.label(spikechan);
eeglabel   = data.label(~spikechan);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for demonstration purpose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = fft_along_rows(x)
y = fft(x, [], 2); % use normal Matlab function to compute the fft along 2nd dimension
