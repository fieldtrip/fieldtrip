function [Sts] = ft_spike_triggeredspectrum(cfg, data, spike)

% FT_SPIKE_TRIGGEREDSPECTRUM computes the Fourier spectrum of the LFP around the spikes.
% Difference to FT_SPIKETRIGGEREDSPECTRUM is that this function takes SPIKE
% input.
%
% Use as  [freq] = ft_spike_triggeredspectrum(cfg,data,spike)
% Important:data.time and spike.trialtime should be
% referenced relative to the same trial trigger times!
%
% Configuration inputs largely match those of ft_freqanalysis_mtmconvol:
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
%     cfg.borderspikes = 'yes' (default) or 'no'. If 'yes', we process the spikes falling at the
%                        border using an LFP that is not centered on the spike.
%
% If the triggered spike leads a spike in another channel, then the angle
% of the Fourier spectrum of that other channel will be negative. 
% % $Id$

% Copyright (C) 2008-2011, Martin Vinck, Robert Oostenveld 
% thanks to Henrique Cabral and Thilo Womelsdorf for testing.

ft_defaults

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();

data = ft_checkdata(data, 'datatype', {'raw'}, 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'foi','t_ftimwin'});
if  isequal(cfg.taper, 'dpss')
  cfg = ft_checkconfig(cfg, 'required', {'tapsmofrq'});
  cfg = ft_checkopt(cfg,'tapsmofrq',{'doublevector', 'doublescalar'});    
end

% get the options
cfg.borderspikes   = ft_getopt(cfg, 'borderspikes','yes');  
cfg.taper          = ft_getopt(cfg, 'taper','hanning');  
cfg.spikechannel   = ft_getopt(cfg,'spikechannel', 'all');
cfg.channel        = ft_getopt(cfg,'channel', 'all');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'taper',{'char', 'function_handle'});
cfg = ft_checkopt(cfg,'borderspikes','char',{'yes', 'no'});
cfg = ft_checkopt(cfg,'t_ftimwin',{'doublevector', 'doublescalar'});
cfg = ft_checkopt(cfg,'foi',{'doublevector', 'doublescalar'});
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'channel', {'cell', 'char', 'double'});

% ensure that all cfg options were valid, and no typos were made

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
nspikesel        = length(cfg.spikechannel); % number of spike channels

% already determine the spectra for the DC, by feeding in zeros in ft_specest_mtmconvol
fsample = 1./ (data.time{1}(2)-data.time{1}(1));
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
timeaxis = -max(numsmp):max(numsmp);
timeaxis = timeaxis./fsample;

% just run once to get the DC spectrum, this creates some additional processing time but quite negligible
keys = ft_cfg2keyval(tmpcfg);
[spec_dc] = ft_specest_mtmconvol(ones(1,length(timeaxis)),timeaxis,keys{:});
spec_dc = nanmean(spec_dc,1); % to average across tapers if multiple tapers were called for

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
  tmpcfg.timeoi    = data.time{iTrial};
  keys = ft_cfg2keyval(tmpcfg);
  [spec,ntapers,freqoi,timeoi] = ft_specest_mtmconvol(data.trial{iTrial},data.time{iTrial},keys{:}); 
  spec = nanmean(spec,1); % for averaging across tapers
  
  % do the DC removal now - this is critical not to introduce peak at low frequencies
  fsel  = find(numsmp<length(data.time{iTrial}));% the frequencies that we can process given the length of the LFP            
  for iFreq = fsel(:)'
      mva = zeros(size(data.trial{iTrial})); % moving average to compute the DC
      for iChan = chansel(:)' 
        % we concove with a kernel that sums to 1 to get the moving average, this is the DC at any time-point
          mva(iChan,:)    = conv(data.trial{iTrial}(iChan,:), ones(1,numsmp(iFreq))./numsmp(iFreq),'same');
      end
      % now simply subtract the DC spectrum from the computed LFP spectrum
      newspec = reshape(spec(1,:,iFreq,:),[nchansel,length(tmpcfg.timeoi)]) - mva.*spec_dc(1,1,iFreq); %subtract DC spec
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
        beg = 1+(numsmp(iFreq)-1)/2; % this is the border
        ed  = length(timeoi) - beg;

        vldSpikes    = unitsmp{iUnit}>=beg & unitsmp{iUnit}<=ed; % determine the valid spikes
        earlySpikes  = unitsmp{iUnit}<beg;
        lateSpikes   = unitsmp{iUnit}>ed;
        if any(earlySpikes)
            rephaseMat(earlySpikes,1,iFreq) = exp(1i*2*pi*freq * (unittime{iUnit}(earlySpikes) - timeoi(beg)) );                    
            spectrum{iTrial}(earlySpikes,:,:) = repmat(spec(:,:,:,beg),[sum(earlySpikes) 1 1]);                    
        end
        if any(lateSpikes) 
            rephaseMat(lateSpikes,1,iFreq)  = exp(1i*2*pi*freq * (unittime{iUnit}(lateSpikes) - timeoi(ed)) );
            spectrum{iTrial}(lateSpikes,:,:) = repmat(spec(:,:,:,ed),[sum(lateSpikes) 1 1]);                
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
Sts.label          = data.label(chansel);
Sts.freq           = freqoi;
Sts.dimord         = 'rpt_chan_freq';
Sts.spikechannel   = spikelabel(spikesel);
for iUnit = 1:nspikesel
    Sts.fourierspctrm{iUnit}  = cat(1, spectrum{iUnit,:});
    Sts.time{iUnit}           = cat(1, spiketime{iUnit,:});  
    Sts.trial{iUnit}          = cat(2, spiketrial{iUnit,:})'; 
end
Sts.trialtime = spike.trialtime;

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, ind] = dbstack;
  cfg.version.name = st(ind);
end
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
Sts.cfg = cfg;

  