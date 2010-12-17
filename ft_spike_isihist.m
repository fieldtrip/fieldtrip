function [isih] = ft_spike_isihist(cfg,spike)

% FT_SPIKE_ISIHSIT computes the interspike interval distribution and statistics.
%
% Use as
%   [ISIH] = ft_spike_isihist(CFG,SPIKE)
%
% Configurations options:
%
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
%   Outputs specific for isihist:
%   isih.avg             = nUnits-by-nBins interspike interval histogram
%   isih.time            = bincenters corresponding to ISIH.ISIH.
%   isih.gammashape      = shape parameter of fitted gamma distribution.
%   isih.gammascale      = scale parameterof fitted gamma distribution.
%   isih.coeffvar        = coefficient of variation (std/mean of isi).
%   isih.isi             = 1-by-nUnits cell with interval to next spike per spike.
%   isih.label

% Copyright (C) 2010, Martin Vinck

if nargin~=2, error('MATLAB:ft_spike_isihist:nargin','Two input arguments required'), end

% check the configuration inputs and enter the defaults
defaults.spikechannel   = {'all'};  
defaults.trials         = {'all'}; 
defaults.bins           = {0:0.002:0.1};
defaults.latency        = {'maxperiod'};              
defaults.outputunit     = {'spikecount' 'proportion'};            
defaults.gammafit       = {'yes' 'no'};
defaults.keeptrials     = {'yes' 'no'};
cfg = ft_spike_sub_defaultcfg(cfg,defaults);

% detect whether the format of spike is correct
hasAllFields = all(isfield(spike, {'time', 'trial', 'trialtime', 'label'}));
if ~hasAllFields, error('MATLAB:ft_spike_isihist:wrongStructInput',...
      'input SPIKE should be struct with .timestamp, .trial, .time, .label fields')
end

% check whether all are of right format
correctInp = iscell(spike.time) & iscell(spike.trial) & iscell(spike.label) & size(spike.trialtime,2)==2;
if ~correctInp, error('MATLAB:ft_spike_isihist:wrongStructInput',...
    '.timestamp, .trial and .label should be cell arrays, time should be nTrials-by-2 matrix')
end

% get the number of trials or change DATA according to cfg.trials
if  strcmp(cfg.trials,'all')    
    cfg.trials = 1:size(spike.trialtime,1);
elseif islogical(cfg.trials)
    cfg.trials = find(cfg.trials); 
elseif ~isrealvec(cfg.trials);
    error('MATLAB:ft_spike_isihist:cfg:trials:wrongInput',...
    'cfg.trials should be logical or numerical selection or string "all"');  
end
cfg.trials = sort(cfg.trials);
nTrials    = length(cfg.trials);

cfg.channel = ft_channelselection(cfg.spikechannel, spike.label);
spikesel    = match_str(spike.label, cfg.channel);
nUnits      = length(spikesel); % number of spike channels
if nUnits==0, error('MATLAB:ft_spike_isihist:cfg:spikechannel:noSpikeChanSelected',...
                    'No spikechannel selected by means of cfg.spikechannel'); 
end

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

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
% remember the configuration details of the input data
try, cfg.previous = spike.cfg; end
% remember the exact configuration details in the output 
isih.cfg     = cfg;

