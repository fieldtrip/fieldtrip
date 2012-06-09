function [par] = ft_spike_isihist(cfg,spike)

% FT_SPIKE_ISI_PAR computes parameters on the inter spike interval
% distribution
%
% The input SPIKE should be organised as the spike or the raw datatype, obtained from
% FT_SPIKE_MAKETRIALS or FT_PREPROCESSING (in that case, conversion is done
% within the function)
%
% Use as
%   [par] = ft_spike_isi_par(cfg, spike)
%
% Configurations:
%   cfg.method  = string, one of
%      'gamfit'      : returns [shape scale] for gamma distribution fit
%      'coeffvar'    : coefficient of variation (sd / mean)
%      'lv'          : Shinomoto's Local Variation measure (2009)
%
%   cfg.spikechannel     = string or index of single spike channel to trigger on (default = 'all')
%                          See FT_CHANNELSELECTION for details
%   cfg.trials           = numeric selection of trials (default = 'all')
%   cfg.latency          = [begin end] in seconds, 'max' (default), 'min', 'prestim'(t<=0), or
%                          'poststim' (t>=0).
%                          If 'max', we use all available latencies.
%                          If 'min', we use only the time window contained by all trials.
%                          If 'prestim' or 'poststim', we use time to or from 0.
% Outputs:
%   isih.avg             = nUnits-by-nBins interspike interval histogram
%   isih.time            = bincenters corresponding to isih.avg
%   isih.isi             = 1-by-nUnits cell with interval to next spike per spike.
%   isih.label           = 1-by-nUnits cell array with labels

% Copyright (C) 2010, Martin Vinck
%
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% check if data is of proper format 
spike = ft_checkdata(spike,'datatype', 'spike', 'feedback', 'yes');

% get the default options
cfg.spikechannel = ft_getopt(cfg,'spikechannel', 'all');
cfg.trials       = ft_getopt(cfg,'trials', 'all');
cfg.latency      = ft_getopt(cfg,'latency','maxperiod');
cfg.method       = ft_getopt(cfg,'method', 'coeffvar');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'latency', {'char', 'ascenddoublebivector'});
cfg = ft_checkopt(cfg,'trials', {'char', 'doublevector', 'logical'}); 
cfg = ft_checkopt(cfg,'method', 'char', {'gamfit', 'coeffvar', 'lv'});

% only allow the specified options
cfg = ft_checkconfig(cfg,'allowed', {'spikechannel', 'latency', 'trials', 'method'});

% get the number of trials or change DATA according to cfg.trials
if  strcmp(cfg.trials,'all')
  cfg.trials = 1:size(spike.trialtime,1);
elseif islogical(cfg.trials)
  cfg.trials = find(cfg.trials);
end
cfg.trials = sort(cfg.trials);
nTrials    = length(cfg.trials);

cfg.channel = ft_channelselection(cfg.spikechannel, spike.label);
spikesel    = match_str(spike.label, cfg.channel);
nUnits      = length(spikesel); % number of spike channels
if nUnits==0, error('No spikechannel selected by means of cfg.spikechannel'); end

% determine the duration of each trial
begTrialLatency = spike.trialtime(cfg.trials,1);
endTrialLatency = spike.trialtime(cfg.trials,2);

% select the latencies
if strcmp(cfg.latency,'minperiod')
  cfg.latency = [max(begTriallatency) min(endTriallatency)];
elseif strcmp(cfg.latency,'maxperiod')
  cfg.latency = [min(begTrialLatency) max(endTrialLatency)];
elseif strcmp(cfg.latency,'prestim')
  cfg.latency = [min(begTrialLatency) 0];
elseif strcmp(cfg.latency,'poststim')
  cfg.latency = [0 max(endTrialLatency)];
end

out = [];
for iUnit = 1:nUnits
  unitIndx = spikesel(iUnit);
  
  % only select the spikes in the right latencies and trials
  spikeTrials    = spike.trial{unitIndx}(:)';
  spikeTimes     = spike.time{unitIndx}(:)';
  
  % select only the spikes in the window and with the selected trials
  spikesInTrials = ismember(spikeTrials, cfg.trials);
  spikesInWin    = spikeTimes>=cfg.latency(1)&spikeTimes<=cfg.latency(2);
  spikeTimes     = spikeTimes(spikesInTrials & spikesInWin);
  spikeTrials    = spikeTrials(spikesInTrials & spikesInWin);
  
  % find the spikes that jumped to the next trial, replace them with NaNs
  trialJump = logical([1 diff(spikeTrials)]);
  isi = [NaN diff(spikeTimes)];
  isi(trialJump) = NaN;

  switch cfg.method
    case 'coeffvar'      
      out(iUnit) = nanstd(isi)./nanmean(isi);
    case 'gamfit'
      data = isi(~isnan(isi));  % remove the nans from isiSpike
      if ~isempty(data)
        [out(iUnit,:)] = mle(data,'distribution', 'gamma'); % fit a gamma distribution
      else
        out(iUnit,:)   = [NaN NaN];
      end
    case 'lv'
      dIsi     = isi(1:end-1) - isi(2:end);
      sumIsi   = isi(1:end-1) + isi(2:end);
      sl       = ~isnan(dIsi) & ~isnan(sumIsi); % remove if one has nan
      df       = sum(sl) - 1;
      out(iUnit) = (3/df)*sum((dIsi(sl)./sumIsi(sl)).^2);
  end
end

% gather the rest of the results
param = cfg.method;
par.dimord       = 'chan_par';
par.(param) = out;
par.label        = spike.label(spikesel);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous spike
ft_postamble history isih


