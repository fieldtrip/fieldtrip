function [isih] = ft_spike_isihist(cfg,spike)

% FT_SPIKE_ISIHIST computes the interspike interval distribution and statistics.
%
% The input SPIKE should be organised as the spike or the raw datatype, obtained from
% FT_SPIKE_MAKETRIALS or FT_PREPROCESSING (in that case, conversion is done
% within the function)
%
% Use as
%   [isih] = ft_spike_isihist(cfg, spike)
%
% Configurations:
%   cfg.outputunit       = 'spikecount' (default) or 'proportion' (sum of all bins = 1).
%   cfg.spikechannel     = string or index of single spike channel to trigger on (default = 'all')
%                          See FT_CHANNELSELECTION for details
%   cfg.trials           = numeric selection of trials (default = 'all')
%   cfg.bins             = vector of isi bins.
%   cfg.latency          = [begin end] in seconds, 'max' (default), 'min', 'prestim'(t<=0), or
%                          'poststim' (t>=0).
%                          If 'max', we use all available latencies.
%                          If 'min', we use only the time window contained by all trials.
%                          If 'prestim' or 'poststim', we use time to or from 0.
%   cfg.gammafit         = 'yes', or 'no' (default), if 'yes', we fit a gamma
%                          distribution to the raw isi times.
%`
% Outputs:
%   isih.avg             = nUnits-by-nBins interspike interval histogram
%   isih.time            = bincenters corresponding to ISIH.ISIH.
%   isih.gammashape      = shape parameter of fitted gamma distribution.
%   isih.gammascale      = scale parameterof fitted gamma distribution.
%   isih.coeffvar        = coefficient of variation (std/mean of isi).
%   isih.isi             = 1-by-nUnits cell with interval to next spike per spike.
%   isih.label

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
cfg.outputunit   = ft_getopt(cfg, 'outputunit','spikecount');
cfg.spikechannel = ft_getopt(cfg, 'spikechannel', 'all');
cfg.trials       = ft_getopt(cfg, 'trials', 'all');
cfg.latency      = ft_getopt(cfg,'latency','maxperiod');
cfg.keeptrials   = ft_getopt(cfg,'keeptrials', 'yes');
cfg.bins         = ft_getopt(cfg,'bins', 0:0.002:1);
cfg.gammafit     = ft_getopt(cfg,'gammafit', 'yes');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'outputunit','char', {'spikecount', 'proportion'});
cfg = ft_checkopt(cfg,'bins', 'doublevector');
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'latency', {'char', 'ascenddoublebivector'});
cfg = ft_checkopt(cfg,'trials', {'char', 'doublevector', 'logical'}); 
cfg = ft_checkopt(cfg,'keeptrials', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'gammafit', 'char', {'yes', 'no'});

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
elseif ~isrealvec(cfg.latency) || length(cfg.latency)~=2
  error('MATLAB:ft_spike_isihist:cfg:latency',....
    'cfg.latency should be "max", "min", "prestim", "poststim" or 1-by-2 vector')
end
if cfg.latency(1)>cfg.latency(2), error('MATLAB:ft_spike_isihist:incorrectLatencyWindow',...
    'cfg.latency should be a vector in ascending order, i.e., cfg.latency(2)>cfg.latency(1)');
end
if (cfg.latency(1) < min(begTrialLatency)), cfg.latency(1) = min(begTrialLatency);
  warning('MATLAB:ft_spike_isihist:incorrectLatencyWindow',...
    'Correcting begin latency of averaging window');
end
if (cfg.latency(2) > max(endTrialLatency)), cfg.latency(2) = max(endTrialLatency);
  warning('MATLAB:ft_spike_isihist:incorrectLatencyWindow',...
    'Correcting end latency of averaging window');
end
doMax = cfg.latency(1)<=min(begTrialLatency) && cfg.latency(2)>=max(endTrialLatency);

% construct the isi bins, we let it run to maximum trial duration by default
bins  = cfg.bins;
nBins = length(bins)-1;

% compute the interspike interval
keepTrials = strcmp(cfg.keeptrials,'yes') | strcmp(cfg.gammafit,'yes');
if keepTrials, isiSpike = cell(1,nUnits); end  % contains the individual histogram spike times
isihist  = zeros(nUnits,nBins+1); % isi histogram

for iUnit = 1:nUnits
  unitIndx = spikesel(iUnit);
  
  % only select the spikes in the right latencies and trials
  spikeTrials    = spike.trial{unitIndx};
  spikeTimes     = spike.time{unitIndx};
  
  % select only the spikes in the window and with the selected trials
  spikesInTrials = ismember(spikeTrials, cfg.trials);
  if ~doMax
    spikesInWin    = spikeTimes>=cfg.latency(1)&spikeTimes<=cfg.latency(2);
  else
    spikesInWin    = 1;
  end
  spikeTimes     = spikeTimes(spikesInTrials & spikesInWin);
  spikeTrials    = spikeTrials(spikesInTrials & spikesInWin);
  
  % find the spikes that jumped to the next trial, replace them with NaNs
  trialJump = logical([1; diff(spikeTrials)]);
  isi = [NaN; diff(spikeTimes)];
  isi(trialJump) = NaN;
  
  % convert and store the isi
  if keepTrials, isiSpike{iUnit} = isi; end
  isihist(iUnit,:)   = histc(isi,bins);
end

isihist(:,end) = []; % the last number is only an equality to a bin edge
if strcmp(cfg.outputunit,'proportion'),
  isihist = isihist./repmat(nansum(isihist,2),1,size(isihist,2));
end

% coefficient of variation
if  strcmp(cfg.gammafit,'yes')
  [shape, scale] = deal(zeros(1,nUnits));
  for iUnit = 1:nUnits
    data = isiSpike{iUnit}(~isnan(isiSpike{iUnit}));  % remove the nans from isiSpike
    if ~isempty(data)
      [par] = mle(data,'distribution', 'gamma'); % fit a gamma distribution
    else
      par   = [NaN NaN];
    end
    shape(iUnit) = par(1); % get the individual parameters and name them
    scale(iUnit) = par(2);
  end
end

% compute the coefficient of variation of the isi
coeffvar = zeros(1,nUnits);
for iUnit = 1:nUnits
  coeffvar(iUnit) = nanstd(isiSpike{iUnit})/nanmean(isiSpike{iUnit});
end

% gather the rest of the results
if strcmp(cfg.keeptrials,'yes')
  isih.isi       = isiSpike;
end
isih.time        = bins(1:end-1);
isih.avg         = isihist;
isih.dimord      = 'chan_bins';
isih.label       = spike.label(spikesel);
if strcmp(cfg.gammafit,'yes')
  isih.gammaShape = shape;
  isih.gammaScale = scale;
end
isih.coeffvar = coeffvar;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous spike
ft_postamble history isih


