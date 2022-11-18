function [sts] = ft_spiketriggeredspectrum_convol(cfg, data, spike)

% FT_SPIKETRIGGEREDSPECTRUM_CONVOL computes the Fourier spectrum (amplitude and
% phase) of the LFP around the spikes using convolution of the complete LFP traces. A
% phase of zero corresponds to the spike being on the peak of the LFP oscillation. A
% phase of 180 degree corresponds to the spike being in the through of the oscillation.
% A phase of 45 degrees corresponds to the spike being just after the peak in the LFP.
%
% The difference to FT_SPIKETRIGGEREDSPECTRUM_FFT is that this function allows for
% multiple frequencies to be processed with different time-windows per frequency, and
% that FT_SPIKETRIGGEREDSPECTRUM_FFT is based on taking the FFT of a limited LFP
% segment around each spike.
%
% Use as
%   [sts] = ft_spiketriggeredspectrum_convol(cfg,data,spike)
% or 
%   [sts] = ft_spiketriggeredspectrum_convol(cfg,data)
% The spike data can either be contained in the data input or in the spike
% input.
%
% The input DATA should be organised as the raw datatype, obtained from
% FT_PREPROCESSING or FT_APPENDSPIKE.
%
% The input SPIKE should be organised as the spike or the raw datatype, obtained from
% FT_SPIKE_MAKETRIALS or FT_PREPROCESSING, in which case the conversion is done
% within this function.
%
% Important is that data.time and spike.trialtime should be referenced
% relative to the same trial trigger times!
%
% Configurations (following largely FT_FREQNALYSIS with method mtmconvol)
%     cfg.tapsmofrq       = vector 1 x numfoi, the amount of spectral smoothing through
%                           multi-tapering. Note that 4 Hz smoothing means
%                           plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
%     cfg.foi             = vector 1 x numfoi, frequencies of interest
%     cfg.taper           = 'dpss', 'hanning' or many others, see WINDOW (default = 'hanning')
%     cfg.t_ftimwin       = vector 1 x numfoi, length of time window (in
%     seconds)
%     cfg.taperopt        =  parameter that goes in WINDOW function (only
%                           applies to windows like KAISER).
%     cfg.spikechannel    = cell-array with selection of channels (default = 'all')
%                           see FT_CHANNELSELECTION for details
%     cfg.channel         = Nx1 cell-array with selection of channels (default = 'all'),
%                           see FT_CHANNELSELECTION for details
%     cfg.borderspikes    = 'yes' (default) or 'no'. If 'yes', we process the spikes
%                           falling at the border using an LFP that is not centered
%                           on the spike. If 'no', we output NaNs for spikes
%                           around which we could not center an LFP segment.
%     cfg.rejectsaturation= 'yes' (default) or 'no'. If 'yes', we set
%                           EEG segments where the maximum or minimum
%                           voltage range is reached
%                           with zero derivative (i.e., saturated signal) to
%                           NaN, effectively setting all spikes phases that
%                           use these parts of the EEG to NaN. An EEG that
%                           saturates always returns the same phase at all
%                           frequencies and should be ignored.
%
% Note: some adjustment of the frequencies can occur as the chosen time-window may not 
% be appropriate for the chosen frequency.
% For example, suppose that cfg.foi = 80, data.fsample = 1000, and
% cfg.t_ftimwin = 0.0625. The DFT frequencies in that case are 
% linspace(0,1000,63) such that cfg.foi --> 80.645. In practice, this error
% can only become large if the number of cycles per frequency is very
% small and the frequency is high. For example, suppose that cfg.foi = 80
% and cfg.t_ftimwin = 0.0125. In that case cfg.foi-->83.33.
% The error is smaller as data.fsample is larger.
%
% Outputs:
%   sts is a spike structure, containing new fields:
%   sts.fourierspctrm = 1 x nUnits cell-array with dimord spike_lfplabel_freq
%   sts.lfplabel      = 1 x nChan cell-array with EEG labels
%   sts.freq          = 1 x nFreq frequencies. Note that per default, not
%                       all frequencies can be used as we compute the DFT
%                       around the spike based on an uneven number of
%                       samples. This introduces a slight adjustment of the
%                       selected frequencies.
%
%   Note: sts.fourierspctrm can contain NaNs, for example if
%   cfg.borderspikes = 'no', or if cfg.rejectsaturation = 'yes', or if the
%   trial length was too short for the window desired.
%
% WHen using multitapering, the phase distortion is corrected for.
%
% The output STS data structure can be input to FT_SPIKETRIGGEREDSPECTRUM_STAT

% Copyright (C) 2008-2012, Martin Vinck
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


% check if the input data is valid for this function
data  = ft_checkdata(data, 'datatype', {'raw'}, 'feedback', 'yes');
if nargin==3
  spike = ft_checkdata(spike, 'datatype', {'spike'}, 'feedback', 'yes');
end

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'foi','t_ftimwin'});

% get the options
cfg.borderspikes   = ft_getopt(cfg, 'borderspikes','yes');
cfg.taper          = ft_getopt(cfg, 'taper','hanning');
cfg.taperopt       = ft_getopt(cfg, 'taperopt',[]);
cfg.spikechannel   = ft_getopt(cfg,'spikechannel', 'all');
cfg.channel        = ft_getopt(cfg,'channel', 'all');
cfg.rejectsaturation  = ft_getopt(cfg,'rejectsaturation', 'yes');

% ensure that the options are valid
cfg = ft_checkopt(cfg, 'taper',{'char', 'function_handle'});
cfg = ft_checkopt(cfg, 'borderspikes','char',{'yes', 'no'});
cfg = ft_checkopt(cfg, 't_ftimwin',{'doublevector', 'doublescalar'});
cfg = ft_checkopt(cfg, 'foi',{'doublevector', 'doublescalar'});
cfg = ft_checkopt(cfg, 'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg, 'channel', {'cell', 'char', 'double'});
cfg = ft_checkopt(cfg, 'taperopt', {'double','empty'});
cfg = ft_checkopt(cfg, 'rejectsaturation','char', {'yes', 'no'});

if  isequal(cfg.taper, 'dpss')
  cfg = ft_checkconfig(cfg, 'required', {'tapsmofrq'});
  cfg = ft_checkopt(cfg,'tapsmofrq',{'doublevector', 'doublescalar'});
end

cfg = ft_checkconfig(cfg, 'allowed', {'taper', 'borderspikes', 't_ftimwin', 'foi', 'spikechannel', 'channel', 'taperopt', 'rejectsaturation','tapsmofrq'});

% length of tapsmofrq, foi and t_ftimwin should all be matched
if isfield(cfg,'tapsmofrq')
  if length(cfg.tapsmofrq) ~= length(cfg.foi) || length(cfg.foi)~=length(cfg.t_ftimwin)
    error('lengths of cfg.tapsmofrq, cfg.foi and cfg_t_ftimwin should be equal and 1 x nFreqs')
  end
end

% get the spikechannels
if nargin==2
  % autodetect the spikechannels and EEG channels
  [spikechannel, eegchannel] = detectspikechan(data);
  if strcmp(cfg.spikechannel, 'all'), 
    cfg.spikechannel = spikechannel; 
  else
    cfg.spikechannel = ft_channelselection(cfg.spikechannel, data.label);  
    if ~all(ismember(cfg.spikechannel,spikechannel)), error('some selected spike channels appear eeg channels'); end        
  end
  if strcmp(cfg.channel,'all')  
    cfg.channel = eegchannel;
  else
    cfg.channel      = ft_channelselection(cfg.channel, data.label);  
    if ~all(ismember(cfg.channel,eegchannel)), warning('some of the selected eeg channels appear spike channels'); end    
  end    
  tmpcfg = [];
  tmpcfg.channel = cfg.spikechannel;
  data_spk = ft_selectdata(tmpcfg, data);
  tmpcfg.channel = cfg.channel;
  data     = ft_selectdata(tmpcfg, data); % leave only LFP
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
if nspikesel==0,   error('no units were selected');   end
if nchansel==0,    error('no channels were selected'); end

% preallocate
nTrials                             = length(data.trial); % number of trials
[spectrum,spiketime, spiketrial]    = deal(cell(nspikesel,nTrials)); % preallocate the outputs
[unitsmp,unittime,unitshift]        = deal(cell(1,nspikesel));
nSpikes = zeros(1,nspikesel);

% construct the frequency axis and restrict to unique frequencies
if ~isfield(data, 'fsample'), data.fsample = 1/mean(diff(data.time{1})); end
if any(cfg.foi > (data.fsample/2)) 
  error('frequencies in cfg.foi are above Nyquist frequency')
end
numsmp  = round(cfg.t_ftimwin .* data.fsample);
numsmp(~mod(numsmp,2)) = numsmp(~mod(numsmp,2))+1; % make sure we always have uneven samples, since we want the spike in the middle
foi = zeros(1,length(cfg.foi));
for iSmp = 1:length(numsmp)
  faxis         = linspace(0,data.fsample,numsmp(iSmp));
  findx         = nearest(faxis,cfg.foi(iSmp));
  [foi(iSmp)]   = deal(faxis(findx)); % this is the actual frequency used, from the DFT formula
end
%
[cfg.foi,B,C] = unique(foi); % take the unique frequencies from this
cfg.t_ftimwin = cfg.t_ftimwin(B);

% compute the minima and maxima of the data, this is done to remove EEG portions where there are potential saturation effects
if strcmp(cfg.rejectsaturation,'yes')
  [minChan,maxChan] = deal([]);
  for iChan = 1:nchansel
    mx = -inf;
    mn = +inf;
    for iTrial = 1:nTrials
       minTrial = nanmin(data.trial{iTrial}(iChan,:));
       maxTrial = nanmax(data.trial{iTrial}(iChan,:));
       if maxTrial>mx, mx = maxTrial; end
       if minTrial<mn, mn = minTrial; end     
    end
    minChan(iChan) = mn;
    maxChan(iChan) = mx;
  end
end

% compute the spectra
ft_progress('init', 'text',     'Please wait...');
for iTrial = 1:nTrials
  
  % select the spike times for a given trial and restrict to those overlapping with the EEG  
  x = data.time{iTrial};      
  timeBins   = [x x(end)+1/data.fsample] - (0.5/data.fsample);      
  for iUnit = 1:nspikesel
    unitindx   = spikesel(iUnit);
    hasTrial   = spike.trial{unitindx} == iTrial; % find the spikes that are in the trial
    ts         = spike.time{unitindx}(hasTrial); % get the spike times for these spikes
    vld        = ts>=timeBins(1) & ts<=timeBins(end); % only select those spikes that fall in the trial window
    ts         = ts(vld); % timestamps for these spikes
    [ignore,I] = histc(ts,timeBins);    
    if ~isempty(ts)
      ts(I==0 | I==length(timeBins)) = [];
    end
    I(I==0 | I==length(timeBins)) = [];
    unitsmp{iUnit}      = I;
    unittime{iUnit}     = ts(:); % this is for storage in the output structure
    smptime = data.time{iTrial}(unitsmp{iUnit}); % times corresponding to samples
    unitshift{iUnit}    = ts(:) - smptime(:); % shift induced by shifting to sample times, important for high-frequency oscillations
    nSpikes(iUnit)      = length(unitsmp{iUnit});
  end
  
  if ~any(nSpikes),continue,end % continue to the next trial if this one does not contain valid spikes
  
  % set the saturated parts of the data to NaNs
  if strcmp(cfg.rejectsaturation,'yes')
    for iChan = 1:nchansel    
      isSaturated = find(diff(data.trial{iTrial}(chansel(iChan),:))==0)+1;
      remove      = data.trial{iTrial}(chansel(iChan),isSaturated)==maxChan(iChan) | data.trial{iTrial}(chansel(iChan),isSaturated)==minChan(iChan);
      remove      = isSaturated(remove);
      data.trial{iTrial}(chansel(iChan),remove) = NaN;
      if ~isempty(remove)
        fprintf('setting %d points from channel %s in trial %d to NaN\n', length(remove), cfg.channel{iChan}, iTrial);
      end
    end
  end
  
  % preallocate
  nFreqs = length(cfg.foi);
  for iUnit = 1:nspikesel
    if nSpikes(iUnit)>0, spectrum{iUnit,iTrial} = zeros(nSpikes(iUnit),nchansel,nFreqs); end
  end
  
  % process the phases for every chan-freq combination separately to save memory
  for iFreq = 1:nFreqs        
    
    % construct the input options for the sub-function phase_est
    tmpcfg               = cfg;
    tmpcfg.foi           = cfg.foi(iFreq);
    try tmpcfg.tapsmofrq = cfg.tapsmofrq(iFreq);end
    tmpcfg.t_ftimwin     = cfg.t_ftimwin(iFreq);
    if tmpcfg.t_ftimwin>=(data.time{iTrial}(end)-data.time{iTrial}(1))
      warning('time window for frequency %.2f Hz too large for trial %d', tmpcfg.foi, iTrial); 
      spectrum{iUnit,iTrial}(:,:,iFreq) = NaN;
      continue;
    end
    
    % just some message
    if nTrials==1
       ft_progress(iFreq/nFreqs, 'Processing frequency %d from %d', iFreq, nFreqs);    
    else
       ft_progress(iTrial/nTrials, 'Processing trial %d from %d', iTrial, nTrials);    
    end      
    % compute the LFP phase at every time-point
    spec = zeros(length(data.time{iTrial}),nchansel);
    for iChan = 1:nchansel
      [spec(:,iChan),foi, numsmp] = phase_est(tmpcfg,data.trial{iTrial}(chansel(iChan),:),data.time{iTrial}, data.fsample);
    end
    
    % select now the proper phases for every unit
    for iUnit = 1:nspikesel
      
      if nSpikes(iUnit)==0, continue,end
      
      % gather the phases at the samples where we had the spikes
      spectrum{iUnit,iTrial}(:,:,iFreq) = spec(unitsmp{iUnit},:);
      
      % preallocate rephasing vector
      rephase   = ones(nSpikes(iUnit),1); 
      
      % do the rephasing and use proper phases for spikes at borders
      if strcmp(cfg.borderspikes,'yes') 
        
        % use an LFP window not centered around the spike, to use spikes around the trial borders as well       
        beg = 1+(numsmp-1)/2; % find the borders empirically
        ed  = length(data.time{iTrial}) - beg + 1;
        vldSpikes    = unitsmp{iUnit}>=beg & unitsmp{iUnit}<=ed; % determine the valid spikes
        earlySpikes  = unitsmp{iUnit}<beg;
        lateSpikes   = unitsmp{iUnit}>ed;
        if any(earlySpikes)
          rephase(earlySpikes) = exp(1i*2*pi*foi* (unittime{iUnit}(earlySpikes) - data.time{iTrial}(beg))' );
          for iChan = 1:nchansel
            spectrum{iUnit,iTrial}(earlySpikes,iChan,iFreq) = spec(beg,iChan);
          end
        end
        if any(lateSpikes)
          rephase(lateSpikes)  = exp(1i*2*pi*foi * (unittime{iUnit}(lateSpikes) - data.time{iTrial}(ed)) );
          for iChan = 1:nchansel
            spectrum{iUnit,iTrial}(lateSpikes,iChan,iFreq) = spec(ed,iChan);
          end
        end
        if any(vldSpikes)
          rephase(vldSpikes)   = exp(1i*2*pi*foi*unitshift{iUnit}(vldSpikes));
        end
        
      else
        
        % in this case the spikes that fall around borders are set to NaN
        if ~isempty(unitshift{iUnit})
          rephase(1:end) = exp(1i*2*pi*foi*unitshift{iUnit});
        end
        
      end
      
      % Rephase the spectrum, this serves two purposes:
      % 1) correct for potential misallignment of spikes to LFP sampling axis
      % 2) for spikes falling at borders, we need to account for the fact that we 
      % used the instantaneous phase at a different moment in time      
      for iChan = 1:nchansel
        spectrum{iUnit,iTrial}(:,iChan,iFreq) = spectrum{iUnit,iTrial}(:,iChan,iFreq).*rephase;
      end
      
      % Store the spiketimes and spiketrials, constructing a type of SPIKE format
      spiketime{iUnit,iTrial}  = unittime{iUnit};
      spiketrial{iUnit,iTrial} = iTrial*ones(1,nSpikes(iUnit));
    end
  end
end
ft_progress('close');
% collect the results
sts.lfplabel          = data.label(chansel);
sts.freq              = cfg.foi;
sts.label             = spike.label(spikesel);
for iUnit = 1:nspikesel
  sts.fourierspctrm{iUnit}  = cat(1, spectrum{iUnit,:}); 
  spectrum(iUnit,:) = {[]}; % free from the memory
  sts.time{iUnit}           = cat(1, spiketime{iUnit,:});
  sts.trial{iUnit}          = cat(2, spiketrial{iUnit,:})';
end
sts.dimord = '{chan}_spike_lfpchan_freq';
sts.trialtime = spike.trialtime;
  
% do the general cleanup and bookkeeping at the end of the function

ft_postamble previous   data spike
ft_postamble provenance sts
ft_postamble history    sts


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spctrm,foi, numsmp] = phase_est(cfg,dat,time,fsample)

% Phase estimation function

% Determine fsample and set total time-length of data
if nargin<4
  fsample = 1./mean(time(2:end)-time(1:end-1)); % round off errors!
end
numsmp  = round(cfg.t_ftimwin .* fsample);
numsmp(~mod(numsmp,2)) = numsmp(~mod(numsmp,2))+1; % make sure we always have uneven samples, since we want the spike in the middle
faxis         = linspace(0,fsample,numsmp);
findx         =  nearest(faxis,cfg.foi);
[cfg.foi,foi] = deal(faxis(findx)); % this is the actual frequency used, from the DFT formula
timwinSamples = numsmp;

% Compute tapers per frequency, multiply with wavelets and compute their fft
switch cfg.taper
  case 'dpss'
    % create a sequence of DPSS tapers
    taper = double_dpss(timwinSamples, timwinSamples .* (cfg.tapsmofrq ./ fsample))';
    taper = taper(1:(end-1), :); % removing the last taper

    % give error/warning about number of tapers
    if isempty(taper)
      error('%.3f Hz: datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',...
        cfg.foi, timwinSamples/fsample,cfg.tapsmofrq,fsample/timwinSamples);
    elseif size(taper,1) == 1
      warning('using only one taper for specified smoothing for %.2f Hz', cfg.foi)
    end

  case 'sine'
    taper = sine_taper(timwinSamples, timwinSamples .* (cfg.tapsmofrq ./ fsample))'; 
    taper = taper(1:(end-1), :);

  otherwise
    % create a single taper according to the standard window specification 
    if isempty(cfg.taperopt)
      taper = window(cfg.taper, timwinSamples)';
    else       
      try
        taper = window(cfg.taper, timwinSamples,cfg.taperopt)';
      catch
        error('taper option was not appropriate for taper');
      end
    end
end

% do some check on the taper size
if size(taper,1)==numsmp, taper = taper'; end

% normalize taper if there's 1 and all are positive
if size(taper,1)==1 && all(taper>=0)
  taper = numsmp*taper./sum(taper); % magnitude of cosine is now returned
end

%%%% fit linear regression for every datapoint: to remove mean and ramp of signal  
sumKern = ones(1,timwinSamples);
avgKern = sumKern./timwinSamples;
xKern   = timwinSamples:-1:1; % because of definition conv2 function
meanX   = mean(xKern);
sumX    = sum(xKern);

% beta1 =  (sum(x.*y) - sum(x)*sum(y)/n) ./ (sum((x-meanx).^2): standard linear regr formula
beta1 = (conv2(dat(:),xKern(:),'same') - sumX.*conv2(dat(:),sumKern(:),'same')/timwinSamples ) ./ sum((xKern-meanX).^2);
beta0 = conv2(dat(:),avgKern(:),'same') - beta1.*meanX; % beta0 = mean(dat) - beta1*mean(x)  

% DFT formula: basefunctions: cos(2*pi*k*[0:numsmp-1]/numsmp) and sin(2*pi*k*[0:numsmp-1]/numsmp)
% center the base functions such that the peak of the cos function is at the center. See:
% f = findx-1; ax = -(numsmp-1)/2:(numsmp-1)/2; y = cos(2*pi.*ax.*f/numsmp);figure, plot(ax,y,'sr'); 
% y2 = cos(2*pi.*linspace(-(numsmp-1)/2,(numsmp-1)/2,1000).*f/numsmp); 
% hold on, plot(linspace(-(numsmp-1)/2,(numsmp-1)/2,1000),y2,'r-')

% correcting for phase rotation (as with multitapering)
nTapers   = size(taper,1);  
indN      = -(numsmp-1)/2:(numsmp-1)/2;
spctrmRot = complex(zeros(1,1));
for iTaper = 1:nTapers
    coswav    =  taper(iTaper,:).*cos(2*pi*(findx-1)*indN/numsmp);
    sinwav    =  taper(iTaper,:).*sin(2*pi*(findx-1)*indN/numsmp);
    wavelet   = complex(coswav(:), sinwav(:));
    cosSignal = cos(2*pi*(findx-1)*indN/numsmp);
    spctrmRot = spctrmRot + sum(wavelet.*cosSignal(:));
end
phaseCor = angle(spctrmRot);

%%%% compute the spectra
nTapers = size(taper,1);  
indN    = -(numsmp-1)/2:(numsmp-1)/2;
spctrm = complex(zeros(length(dat),1));
for iTaper = 1:nTapers
    coswav  = taper(iTaper,:).*cos(2*pi*(findx-1)*indN/numsmp);
    sinwav  = taper(iTaper,:).*sin(2*pi*(findx-1)*indN/numsmp);
    wavelet = complex(coswav(:), sinwav(:));       
    fftRamp = sum(xKern.*coswav) + 1i*sum(xKern.*sinwav); % fft of ramp with dx/ds = 1 * taper 
    fftDC   = sum(ones(1,timwinSamples).*coswav) + 1i*sum(ones(1,timwinSamples).*sinwav); % fft of unit direct current * taper
    spctrm  = spctrm + (conv_fftbased(dat(:),wavelet) - (beta0*fftDC + beta1.*fftRamp))/(numsmp/2);           
                       % fft                       % mean            %linear ramp      % make magnitude invariant to window length                             
end
spctrm = spctrm./nTapers; % normalize by number of tapers
spctrm = spctrm.*exp(-1i*phaseCor);

% set the part of the spectrum without a valid phase to NaNs
n = (timwinSamples-1)/2;
spctrm(1:n) = NaN; 
spctrm(end-n+1:end) = NaN;

function [taper] = double_dpss(a, b, varargin)
taper = dpss(double(a), double(b), varargin{:});


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

% CONVOLUTION: FFT BASED IMPLEMENTATION
function c = conv_fftbased(a, b)

P = numel(a);
Q = numel(b);
L = P + Q - 1;
K = 2^nextpow2(L);

c = ifft(fft(a, K) .* fft(b, K));
c = c(1:L);
toRm1 = [1:(Q-1)/2];
toRm2 = [(1 + length(c) - (Q-1)/2) : length(c)];
toRm  = [toRm1(:); toRm2(:)];
c(toRm) = [];
