function [rate] = ft_spike_rate(cfg,spike)

% FT_SPIKE_RATE computes the firing rate of spiketrains and their variance
%
% Use as
%   [RATE] = FT_SPIKE_RATE(CFG,SPIKE)
%   
% The input SPIKE should be organised as:   
%   The SPIKE datatype, obtained from FT_SPIKE_DATA2SPIKE or
%   FT_SPIKE_MAKETRIALS.
%
% Configurations options:
%
%   cfg.outputunit       =  - 'rate' (default) If 'rate', we convert the 
%                              output per trial to firing rates (spikes/sec).
%                           - 'spike': we count the number spikes per trial.
%   cfg.spikechannel     =  See FT_CHANNELSELECTION for details
%   cfg.trials           =  - vector of indices (e.g., 1:2:10)
%                           - logical selection of trials (e.g., [1010101010])
%                           - 'all' (default), selects all trials%   cfg.trials                                    
%   cfg.vartriallen      =  'yes' (default) or 'no'. 
%                           If 'yes' - accept variable trial lengths and use all available trials
%                           and the samples in every trial. 
%                           If 'no'  - only select those trials that fully cover the window as
%                           specified by cfg.latency and discard those trials that do not.        
%   cfg.latency          =  - [begin end] in seconds
%                           - 'maxperiod' (default)
%                           - 'minperiod', i.e., the minimal period all trials share
%                           - 'prestim' (all t<=0)
%                           - 'poststim' (all t>=0). 
%   cfg.keeptrials       = 'yes' or 'no' (default). 
%
% The outputs from spike are the following:
%       - rate.trial:        nTrials x nUnits matrix containing the firing rate per unit
%                            and trial
%       - rate.avg:          nTrials array containing the average firing rate per unit
%       - rate.var:          nTrials array containing the variance of firing rates per unit
%       - rate.dof:          nTrials array containing the degree of freedom per unit
%       - rate.label:        nUnits cell array containing the labels of the neuronal units%

% Copyright (C) 2010, Martin Vinck
% removed covariance / correlation here since one may get it very easily with COV or
% corrcoef
% ADD: selection on condition
% ADD: population script, for correlation - principal components - population vector

if nargin~=2, error('ft:spike_rate:nargin','Two input arguments required'), end

% defaults
defaults.outputunit         = {'rate' 'spikecount'};             
defaults.spikechannel       = {'all'}; 
defaults.trials             = {'all'};             
defaults.latency            = {'maxperiod'};              
defaults.vartriallen        = {'yes' 'no'};              
defaults.keeptrials         = {'yes' 'no'};               
cfg = ft_spike_sub_defaultcfg(cfg,defaults);

% detect the format of spike, replace with CHECKDATA eventually, but that needs modification
hasAllFields = all(isfield(spike, {'time', 'trial', 'trialtime', 'label'}));
if ~hasAllFields, error('ft:spike_rate:wrongStructInput',...
      'input spike should be struct with .time, .trial, .trialtime, .label fields')
end

% check whether all are of right format
correctInp = iscell(spike.time) && iscell(spike.trial) && iscell(spike.label) && isrealmat(spike.trialtime) ...
  && size(spike.trialtime,2)==2;
if ~correctInp, error('ft:spike_rate:wrongStructInput',...
    '.time, .trial and .label should be cell arrays, trialtime should be nTrials-by-2 matrix')
end

% get the spikechannels
cfg.channel = ft_channelselection(cfg.spikechannel, spike.label);
spikesel    = match_str(spike.label, cfg.channel);
nUnits      = length(spikesel); % number of spike channels
if nUnits==0, error('ft:spike_rate:cfg:spikechannel:noSpikeChanSelected',...
                    'No spikechannel selected by means of cfg.spikechannel'); 
end

% get the number of trials 
cfg        = trialselection(cfg,spike);
nTrials    = length(cfg.trials); % actual number of trials we use

% determine the duration of each trial
begTrialLatency = spike.trialtime(cfg.trials,1); 
endTrialLatency = spike.trialtime(cfg.trials,2);  

% select the latencies
cfg = latencyselection(cfg,begTrialLatency,endTrialLatency);

% check which trials will be used based on the latency
overlaps      = endTrialLatency>(cfg.latency(1)) & begTrialLatency<(cfg.latency(2));
hasWindow     = ones(nTrials,1);
if strcmp(cfg.vartriallen,'no') % only select trials that fully cover our latency window
    startsLater    = begTrialLatency>cfg.latency(1);
    endsEarlier    = endTrialLatency<cfg.latency(2);
    hasWindow      = ~(startsLater | endsEarlier); % it should not start later or end earlier
    trialDur       = ones(nTrials,1)*(cfg.latency(2)-cfg.latency(1)); 
elseif strcmp(cfg.vartriallen,'yes')
    winBeg      = max([begTrialLatency(:) cfg.latency(1)*ones(nTrials,1)],[],2);
    winEnd      = min([endTrialLatency(:) cfg.latency(2)*ones(nTrials,1)],[],2);
    trialDur    = winEnd-winBeg;
end
cfg.trials         = cfg.trials(overlaps(:) & hasWindow(:));
trialDur           = trialDur(overlaps(:) & hasWindow(:)); % select the durations, we will need this later
nTrials            = length(cfg.trials);
% issue an explicit error if nothing was selected, this is in general indicative of bug
if isempty(cfg.trials), warning('ft:spike_rate:cfg:trials:noneSelected',...
   'No trials were selected in the end, please give us something to analyse'); 
end

% preallocate before computing the psth
keepTrials = strcmp(cfg.keeptrials,'yes');
if  keepTrials, singleTrials = NaN(nTrials,nUnits); end % preallocate single trials with NaNs
[s,ss]   = deal(zeros(nUnits,1));
dof      = nTrials*ones(nUnits,1); % compute degrees of freedom
for iUnit = 1:nUnits 
    unitIndx   = spikesel(iUnit);        
    ts         = spike.time{unitIndx}(:); % get the times
    latencySel = ts>=cfg.latency(1) & ts<=cfg.latency(2); % select timestamps within latency
    trialNums  = spike.trial{unitIndx}(latencySel);% trial indices        

    % use the fact that trial numbers are integers >=1 apart, so we can use histc
    trialBins   = sort([cfg.trials-0.5; cfg.trials+0.5]); 
    trialRate   = histc(trialNums(:),trialBins); 
    trialRate   = trialRate(1:2:end-1); % the uneven bins correspond to the trial integers
    if isempty(trialRate), trialRate = zeros(nTrials,1); end
    
    % convert to firing rates if requested 
    if strcmp(cfg.outputunit,'rate') ,trialRate = trialRate(:)./trialDur(:); end 
    
    % store and compute sum, just fill the single trials up with nans
    if keepTrials, singleTrials(:,iUnit) = trialRate(:); end
    s(iUnit)  = sum(trialRate);
    ss(iUnit) = sum(trialRate.^2);
end

% compute the average rate
rate.avg       = s ./ dof;
rate.var       = (ss - s.^2./dof)./(dof); % since sumrate.^2 ./ dof = dof .* (sumrate/dof).^2
rate.var(dof<2)  = 0;

% gather the rest of the results
rate.dof       = dof;
rate.label     = spike.label(spikesel);   
rate.dimord    = 'chan';
if (strcmp(cfg.keeptrials,'yes'))
  rate.trial = singleTrials;
  rate.dimord = 'rpt_chan';
end
  
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
rate.cfg     = cfg;


%%%%%%%%% SUB FUNCTIONS %%%%%%%%%
function [cfg] = latencyselection(cfg,begTrialLatency,endTrialLatency)

if strcmp(cfg.latency,'minperiod')
      cfg.latency = [max(begTrialLatency) min(endTrialLatency)];
elseif strcmp(cfg.latency,'maxperiod')
      cfg.latency = [min(begTrialLatency) max(endTrialLatency)];
elseif strcmp(cfg.latency,'prestim')
      cfg.latency = [min(begTrialLatency) 0];
elseif strcmp(cfg.latency,'poststim')
      cfg.latency = [0 max(endTrialLatency)];
elseif ~isrealvec(cfg.latency)||length(cfg.latency)~=2
  error('ft:spike_rate:cfg:latency',...
    'cfg.latency should be "max", "min", "prestim", "poststim" or 1-by-2 numerical vector');
end
if cfg.latency(1)>cfg.latency(2), error('ft:spike_rate:incorrectLatencyWindow',...
    'cfg.latency should be a vector in ascending order, i.e., cfg.latency(2)>cfg.latency(1)');
end

% check whether the time window fits with the data
if (cfg.latency(1) < min(begTrialLatency)), cfg.latency(1) = min(begTrialLatency); 
  warning('ft:spike_rate:incorrectLatencyWindow','%s %2.2f',...
          'Correcting begin latency of averaging window to ', cfg.latency(1));
end
if (cfg.latency(2) > max(endTrialLatency)), cfg.latency(2) = max(endTrialLatency);
  warning('ft:spike_rate:incorrectLatencyWindow','%s %2.2f',...
          'Correcting end latency of averaging window to ',cfg.latency(2));
end

%%%

function [cfg] = trialselection(cfg,spike)

% get the number of trials or change DATA according to cfg.trials
nTrials = size(spike.trialtime,1);     
if  strcmp(cfg.trials,'all')    
    cfg.trials = 1:nTrials;
elseif islogical(cfg.trials)
    cfg.trials = find(cfg.trials); 
elseif ~isrealvec(cfg.trials);
    error('ft:spike_rate:cfg:trials:unknownOption',...
    'cfg.trials should be logical or numerical selection or string "all"');  
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>nTrials, error('ft:spike_rate:cfg:trials:maxExceeded',...
   'maximum trial number in cfg.trials should not exceed length of DATA.trial')
end
if isempty(cfg.trials), error('ft:spike_rate:cfg:trials:noneSelected',...
   'No trials were selected by you, rien ne va plus'); 
end


