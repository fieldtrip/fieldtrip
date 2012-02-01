function [Sts] = ft_spike_triggeredspectrum(cfg, data, spike)

% FT_SPIKE_TRIGGEREDSPECTRUM computes the Fourier spectrum of the LFP
% around the spikes. The difference to FT_SPIKETRIGGEREDSPECTRUM is that
% this function allows for multiple frequencies to be processed with
% different time-windows per frequency.
%
% The input SPIKE should be organised as the spike or the raw datatype, obtained from
% FT_SPIKE_MAKETRIALS or FT_PREPROCESSING (in that case, conversion is done
% within the function)
%
% The input DATA should be organised as the raw datatype, obtained from
% FT_PREPROCESSING
%
% Use as
%   [freq] = ft_spike_triggeredspectrum(cfg,data,spike)
%
% Important is that data.time and spike.trialtime should be referenced
% relative to the same trial trigger times!
%
% Configurations:
%     cfg.tapsmofrq    = vector 1 x numfoi, the amount of spectral smoothing through
%                        multi-tapering. Note that 4 Hz smoothing means
%                        plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
%     cfg.foi          = vector 1 x numfoi, frequencies of interest
%     cfg.taper        = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%     cfg.t_ftimwin    = vector 1 x numfoi, length of time window (in seconds)
%     cfg.spikechannel = cell-array with selection of channels (default = 'all')
%                        see FT_CHANNELSELECTION for details
%     cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                        see FT_CHANNELSELECTION for details
%     cfg.borderspikes = 'yes' (default) or 'no'. If 'yes', we process the spikes
%                        falling at the border using an LFP that is not centered
%                        on the spike.
%     cfg.taperopt     = parameter that goes in WINDOW function (only
%                        applies to windows like KAISER)
% If the triggered spike leads a spike in another channel, then the angle
% of the Fourier spectrum of that other channel will be negative.

% Copyright (C) 2008-2011, Martin Vinck, Robert Oostenveld
% thanks to Henrique Cabral and Thilo Womelsdorf for testing.
%
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% check if the input data is valid for this function
data  = ft_checkdata(data, 'datatype', {'raw'}, 'feedback', 'yes');
spike = ft_checkdata(spike, 'datatype', {'spike'}, 'feedback', 'yes');

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'foi','t_ftimwin'});
if  isequal(cfg.taper, 'dpss')
  cfg = ft_checkconfig(cfg, 'required', {'tapsmofrq'});
  cfg = ft_checkopt(cfg,'tapsmofrq',{'doublevector', 'doublescalar'});
end

% get the options
cfg.borderspikes   = ft_getopt(cfg, 'borderspikes','yes');
cfg.taper          = ft_getopt(cfg, 'taper','hanning');
cfg.taperopt       = ft_getopt(cfg, 'taperopt',[]);
cfg.spikechannel   = ft_getopt(cfg,'spikechannel', 'all');
cfg.channel        = ft_getopt(cfg,'channel', 'all');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'taper',{'char', 'function_handle'});
cfg = ft_checkopt(cfg,'borderspikes','char',{'yes', 'no'});
cfg = ft_checkopt(cfg,'t_ftimwin',{'doublevector', 'doublescalar'});
cfg = ft_checkopt(cfg,'foi',{'doublevector', 'doublescalar'});
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'channel', {'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'taperopt', {'double','empty'});
  
% length of tapsmofrq, foi and t_ftimwin should all be matched
if isfield(cfg,'tapsmofrq')
  if length(cfg.tapsmofrq) ~= length(cfg.foi)
    error('length of cfg.tapsmofrq should equal length of cfg.foi')
  end
end

if length(cfg.foi)~=length(cfg.t_ftimwin),
  error('length of cfg.t_ftimwin should equal length of cfg.foi')
end

% determine the channels to be averaged
cfg.channel = ft_channelselection(cfg.channel, data.label);
chansel     = match_str(data.label, cfg.channel); % selected channels
nchansel    = length(cfg.channel);                % number of channels

% get the spikechannels
spikelabel = spike.label;
cfg.spikechannel = ft_channelselection(cfg.spikechannel, spikelabel);
spikesel         = match_str(spikelabel, cfg.spikechannel);
nspikesel        = length(spikesel); % number of spike channels

% get the options that will go in the ft_specest_mtmconvol function
fsample = 1./ mean(diff(data.time{1})); % assumes constancy
tmpcfg = [];
tmpcfg.taper  = cfg.taper;
tmpcfg.freqoi = cfg.foi;
tmpcfg.timwin = cfg.t_ftimwin;
if strcmp(cfg.taper,'dpss'), tmpcfg.tapsmofrq = cfg.tapsmofrq; end
tmpcfg.verbose = 0;

% compute the number of samples, and construct a timeaxis for the DC spectrum
numsmp = round(tmpcfg.timwin .* fsample);
numsmp(~mod(numsmp,2)) = numsmp(~mod(numsmp,2))+1; % make sure we always have uneven samples, since we want the spike in the middle
tmpcfg.timwin = numsmp ./ fsample;
tmpcfg.timeoi    = 0;
tmpcfg.polyorder = -1; % since we do this ourselves in a time-resolved way (it has to)
tmpcfg.taperopt  = cfg.taperopt;

% preallocate
nTrials    = length(data.trial); % number of trials
[spectrum,spiketime, spiketrial]    = deal(cell(nspikesel,nTrials)); % preallocate the outputs
[unitsmp,unittime,unitshift]        = deal(cell(1,nspikesel));
nSpikes = zeros(1,nspikesel);

% compute the spectra
for iTrial = 1:nTrials
  for iUnit = 1:nspikesel
    unitindx = spikesel(iUnit);
    hasTrial = spike.trial{unitindx} == iTrial; % find the spikes that are in the trial
    ts       = spike.time{unitindx}(hasTrial); % get the spike times for these spikes
    vld      = ts>=data.time{iTrial}(1) & ts<=data.time{iTrial}(end); % only select those spikes that fall in the trial window
    ts       = ts(vld); % timestamps for these spikes\
    unitsmp{iUnit} = zeros(1,length(ts));
    for iSpike  = 1:length(ts) % changed this to cope with an sampling frequency that is not a constant but with some jitter
      unitsmp{iUnit}(iSpike) = nearest(data.time{iTrial}, ts(iSpike));
    end
    unittime{iUnit}   = ts(:); % this is for storage in the output structure
    smptime = data.time{iTrial}(unitsmp{iUnit}); % times corresponding to samples
    unitshift{iUnit}    = ts(:) - smptime(:); % shift induced by shifting to sample times, important for high-frequency oscillations
    nSpikes(iUnit)   = length(unitsmp{iUnit});
  end
  if ~any(nSpikes),continue,end % continue to the next trial if this one does not contain valid spikes
  
  % compute the lfp spectrum for all the timepoints and do the selection afterwards
  tmpcfg.timeoi    = 'all';
  keys = ft_cfg2keyval(tmpcfg);
  
  % to compute the DC which must be subtracted
  len = size(data.trial{iTrial},2);
  [specDC,ntapers,freqoi,timeoi] = ft_specest_mtmconvol(ones(1,len),data.time{iTrial},keys{:});
  specDC = nanmean(specDC,1);
  
  % to compute the actual spectra
  [spec,ntapers,freqoi,timeoi] = ft_specest_mtmconvol(data.trial{iTrial}(chansel,:),data.time{iTrial},keys{:});
  spec = nanmean(spec,1); % for averaging across tapers
  
  % do the DC removal now - this is critical not to introduce peak at low frequencies
  fsel  = find(numsmp<length(data.time{iTrial}));% the frequencies that we can process given the length of the LFP
  for iFreq = fsel(:)'
    mva = zeros(nchansel,size(data.trial{iTrial},2)); % moving average to compute the DC
    for iChan = chansel(:)'
      % we concove with a kernel that sums to 1 to get the moving average, this is the DC at any time-point
      mva(iChan,:)    = conv2(data.trial{iTrial}(iChan,:)', ones(1,numsmp(iFreq))'./numsmp(iFreq),'same');
    end
    % now simply subtract the DC spectrum from the computed LFP spectrum
    newspec = reshape(spec(1,:,iFreq,:),[nchansel,length(timeoi)]) - mva.*repmat(squeeze(specDC(1,1,iFreq,:))',[nchansel,1]); %subtract DC spec
    spec(1,:,iFreq,:) = reshape(newspec,[1 nchansel 1 length(timeoi)]);
  end
    
  nFreqs = length(freqoi);
  for iUnit = 1:nspikesel
    fprintf('processing trial %d of %d (%d spikes) for unit %s \n', iTrial, nTrials, nSpikes(iUnit), spikelabel{iUnit});
    spectrum{iUnit,iTrial} = permute(shiftdim(spec(:,:,:,unitsmp{iUnit}),1),[3 1 2]); %%changed by Henrique (before:permute(shiftdim(spec(:,:,:,unitsmp{iUnit}),1),[3 1 2])
    rephaseMat   = ones(nSpikes(iUnit),1,nFreqs);
    for iFreq = fsel(:)'
      freq = freqoi(iFreq);
      
      if strcmp(cfg.borderspikes,'yes') % if yes, we use an LFP window not centered around the spike
        
        beg = 1+(numsmp(iFreq)-1)/2; % find the borders empirically
        ed  = length(timeoi) - beg;
        
        vldSpikes    = unitsmp{iUnit}>=beg & unitsmp{iUnit}<=ed; % determine the valid spikes
        earlySpikes  = unitsmp{iUnit}<beg;
        lateSpikes   = unitsmp{iUnit}>ed;
        if any(earlySpikes)
          rephaseMat(earlySpikes,1,iFreq) = exp(1i*2*pi*freq * (unittime{iUnit}(earlySpikes) - timeoi(beg)) );
          spectrum{iUnit,iTrial}(earlySpikes,:,iFreq) = repmat(spec(:,:,iFreq,beg),[sum(earlySpikes) 1 1]);
        end
        if any(lateSpikes)
          rephaseMat(lateSpikes,1,iFreq)  = exp(1i*2*pi*freq * (unittime{iUnit}(lateSpikes) - timeoi(ed)) );
          spectrum{iUnit,iTrial}(lateSpikes,:,iFreq) = repmat(spec(:,:,iFreq,ed),[sum(lateSpikes) 1 1]);
        end
        if any(vldSpikes)
          rephaseMat(vldSpikes,1,iFreq)   = exp(1i*2*pi*freq*unitshift{iUnit}(vldSpikes));
        end
      else
        if ~isempty(unitshift{iUnit})
          rephaseMat(1:end,1,iFreq) = exp(1i*2*pi*freq*unitshift{iUnit});
        end
      end
    end
    rephaseMat = repmat(rephaseMat,[1 nchansel 1]);
    spectrum{iUnit,iTrial} = spectrum{iUnit,iTrial}.*rephaseMat;
    
    % store the spiketimes and spiketrials, constructing a type of SPIKE format
    spiketime{iUnit,iTrial}  = unittime{iUnit};
    spiketrial{iUnit,iTrial} = iTrial*ones(1,nSpikes(iUnit));
  end
end

% collect the results
Sts.lfplabel          = data.label(chansel);
Sts.freq           = freqoi;
Sts.dimord         = 'rpt_chan_freq';
Sts.label   = spikelabel(spikesel);
for iUnit = 1:nspikesel
  Sts.fourierspctrm{iUnit}  = cat(1, spectrum{iUnit,:});
  Sts.time{iUnit}           = cat(1, spiketime{iUnit,:});
  Sts.trial{iUnit}          = cat(2, spiketrial{iUnit,:})';
end
Sts.fourierspctrmdimord = '{chan}_spike_lfpchan_freq';
Sts.trialtime = spike.trialtime;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous data spike
ft_postamble history Sts

