function [psth] = ft_spike_psth(cfg,spike)
% FT_SPIKE_PSTH computes the peristimulus histogram of spiketrains.
%
% Use as
%   [PSTH] = FT_SPIKE_PSTH(CFG,SPIKE)
%
% The input SPIKE should be organised as:   
%   The SPIKE datatype, obtained from FT_SPIKE_DATA2SPIKE or
%   FT_SPIKE_MAKETRIALS.
%
% Configurations options:
%
%   cfg.binsize          =  [binsize] in sec (default = 0.025 sec); 
%   cfg.outputunit       =  - 'rate' (default) If 'rate', we convert the 
%                              output per bin to firing rates (spikes/sec).
%                           - 'spike': we count the number spikes per bin.
%   cfg.spikechannel     =  See FT_CHANNELSELECTION for details.
%   cfg.trials           =  - vector of indices (e.g., 1:2:10)
%                           - logical selection of trials (e.g., [1010101010])
%                           - 'all' (default), selects all trials
%   cfg.vartriallen      =  - 'yes' (default)
%                             Accept variable trial lengths and use all available
%                             trials and the samples in every trial. Missing values will be
%                             ignored in the computation of the average and the variance and
%                             stored as NaNs in the output PSTH.TRIAL. 
%                           - 'no'  
%                             Only select those trials that fully cover the window as specified
%                             by cfg.latency and discard those trials that do not.
%   cfg.latency          =  - [begin end] in seconds
%                           - 'maxperiod' (default), i.e., maximum period available
%                           - 'minperiod', i.e., the minimal period all trials share
%                           - 'prestim' (all t<=0)
%                           - 'poststim' (all t>=0). 
%   cfg.keeptrials       =  'yes' or 'no' (default).
%
%   Outputs:
%     Psth is a structure similar to FT_TIMELOCKANALYSIS or FT_SPIKE_DENSITY.
%     Psth.time        = center histogram bin points
%	    Psth.fsample 		 = 1/binsize; 
%     Psth.avg         = contains average PSTH per unit
%     Psth.trial       = contains PSTH per unit per trial 
%     Psth.var         = contains variance of PSTH per unit across trials
%   
%     PSTH can be fed into FT_TIMELOCKSTATISTICS, into FT_SPIKE_PLOT_PSTH (for
%     plotting the PSTH), into FT_SPIKE_JPSTH (calculating the joint psth), into
%     FT_SPIKE_PLOT_RASTER as cfg.topdata, in which we plot it above a rasterplot.
%

%  Copyright (C) 2010, Martin Vinck
%  TO ADD: binsize should be optimally determined based on
%  standard techniques to compute binsize based on bandwidth. Same for FT_SPIKE_DENSITY

if nargin~=2, error('ft:spike_psth:nargin','Two input arguments required'), end
   
% check configuration inputs and enter the defaults, should be standardized within fieldtrip
defaults.binsize           = {0.025};      % sec              
defaults.outputunit        = {'rate' 'spikecount'};             
defaults.spikechannel      = {'all'};              
defaults.trials            = {'all'};              
defaults.latency           = {'maxperiod'};              
defaults.vartriallen       = {'yes' 'no'};              
defaults.keeptrials        = {'yes' 'no'};               
cfg = ft_spike_sub_defaultcfg(cfg,defaults);

% detect whether the format of spike is correct, 
% should - when the spike format solidifies, in the end be done with checkdata in fieldtrip 
hasAllFields = all(isfield(spike, {'trial', 'time', 'trialtime' 'label'}));
if ~hasAllFields, error('ft:spike_psth:wrongStructInput',...
      'input SPIKE should be structure  with .trial, .time, .trialtime .label fields')
end

% check whether all are of right format
correctInp = iscell(spike.time) && iscell(spike.trial) && iscell(spike.label) && isrealmat(spike.trialtime) ...
  && size(spike.trialtime,2)==2;
if ~correctInp, error('ft:spike_psth:wrongStructInput',...
    '.time, .trial and .label should be cell arrays, .trialtime should be nTrials-by-2 matrix')
end

% get the number of trials or change DATA according to cfg.trials
cfg        = trialselection(cfg,spike);
 
% select the unit - this should be done with channelselection function 
cfg.channel = ft_channelselection(cfg.spikechannel, spike.label);
spikesel    = match_str(spike.label, cfg.channel);
nUnits      = length(spikesel);
if nUnits==0, error('ft:spike_psth:cfg:spikechannel:noSpikeChanSelected',...
                    'No spikechannel selected by means of cfg.spikechannel'); 
end
  
% determine the duration of each trial - we assume N by 2, see error check before
begTrialLatency = spike.trialtime(cfg.trials,1); % remember: already selected on trial here
endTrialLatency = spike.trialtime(cfg.trials,2);  
trialDur 		= endTrialLatency - begTrialLatency;
% while we could simply deselect trials with trialtime field messed up, this may detect bugs
if any(trialDur<0), error('MATLAB:spike_psth:SPIKE.time:reversedOrder',...
  'SPIKE.time(:,1) should preceed SPIKE.time(:,2), your SPIKE.time field is messed up');
end
 
% select the latencies, use the same modular function in all the scripts
cfg = latencyselection(cfg,begTrialLatency,endTrialLatency);

% do some error checking on the binsize
if cfg.binsize<=0 || cfg.binsize>(cfg.latency(2)-cfg.latency(1)),
   error('ft:spike_psth:cfg:binsize:wrongInput',...
  'cfg.binsize should be greater than zero and not exceed the trialduration');
end

% end of trial should be late enough, beginning should be early enough
fullDur       = trialDur>=cfg.binsize; % trials which have full duration
overlaps      = endTrialLatency>=(cfg.latency(1)+cfg.binsize) & begTrialLatency<=(cfg.latency(2)-cfg.binsize);
if strcmp(cfg.vartriallen,'no') % only select trials that fully cover our latency window
    startsLater    = single(begTrialLatency)>single(cfg.latency(1));
    endsEarlier    = single(endTrialLatency)<single(cfg.latency(2));
    hasWindow      = ~(startsLater | endsEarlier); % it should not start later or end earlier
else
    hasWindow     = ones(1,length(cfg.trials));
end
trialSel           = fullDur(:) & overlaps(:) & hasWindow(:);
cfg.trials         = cfg.trials(trialSel); % note that endTrialLatency was of length cfg.trials 
if isempty(cfg.trials), % it should be explicitly tested that selecting no trial gives no error here
  warning('ft:spike_psth:cfg:trials','No trials were selected after latency selection'); 
end
nTrials         = length(cfg.trials);
begTrialLatency = begTrialLatency(trialSel); % note that begTrialLatency was of length cfg.trials here
endTrialLatency = endTrialLatency(trialSel);

% create the bins, this might harbour some inaccuracy if the latency window is not an
% integer multiple; we start arbitrarily at cfg.latency(1)
bins  = cfg.latency(1) : cfg.binsize : cfg.latency(end);
nBins = length(bins);

% preallocate before computing the psth, otherwise, compute stuff 'on the run'
if strcmp(cfg.keeptrials,'yes'), singleTrials = NaN(nTrials,nUnits,nBins-1); end 
[s,ss]   = deal(zeros(nUnits, nBins-1)); %sum, and sum of squared in the loop, to get mean & var

% preallocate and compute degrees of freedom
allStartEarlier =  all(begTrialLatency<=cfg.latency(1));
allEndLater     =  all(endTrialLatency>=cfg.latency(2));
if  ~ (allStartEarlier && allEndLater), 
  dof = zeros(1,nBins-1); 
else
  dof = ones(1,nBins-1)*nTrials; % if all exceed latency, then this is dof  
end

% compute the peristimulus histogram, different algorithms per data type
for iTrial = 1:nTrials % nTrials redefined on line 149, only the selected trials
    origTrial = cfg.trials(iTrial);     
    if  ~ (allStartEarlier && allEndLater) % select bins and count dof + 1
        binSel = begTrialLatency(iTrial)<=bins(1:end-1) & endTrialLatency(iTrial)>=bins(2:end);
        dof(binSel)      = dof(binSel) + 1;
    else
        binSel           = 1:(nBins-1); % always deselect the last bin
    end
 
    for iUnit = 1:nUnits                        
        unitIndx      = spikesel(iUnit); % select the unit                      
        spikesInTrial = (spike.trial{unitIndx}==origTrial); % get spikes in trial
        spikeTimes    = spike.time{unitIndx}(spikesInTrial);
     
        % compute the psth 
        trialPsth   = histc(spikeTimes(:), bins); % we deselect the last bin per default
        trialPsth   = trialPsth(:)'; % force into row vec
        if isempty(trialPsth), trialPsth = zeros(1,length(bins)); end
        
        % convert to firing rates if requested, with spikecount do nothing
        if strcmp(cfg.outputunit,'rate'), trialPsth = trialPsth/cfg.binsize; end

        % compute the sum and the sum of squares for the var and the mean on the fly
        s(iUnit,binSel)  = s(iUnit,binSel)   + trialPsth(binSel); % last bin is single value
        ss(iUnit,binSel) = ss(iUnit,binSel)  + trialPsth(binSel).^2;

        if strcmp(cfg.keeptrials,'yes'), 
          singleTrials(iTrial,iUnit,binSel) = trialPsth(binSel);
        end                
    end
end

% get the results
dof            = dof(ones(nUnits,1),:);
psth.avg       = s ./ dof;
psth.var       = (ss - s.^2./dof)./(dof-1); % since sumPsth.^2 ./ dof = dof .* (sumPsth/dof).^2
psth.dof       = dof; % combined with psth.var we can get SEM
psth.fsample   = 1/(cfg.binsize);   % might be more compatible with feeding psth in other funcs
psth.time      = bins(1:end-1) + 0.5*cfg.binsize; % use .time name
psth.label     = spike.label(spikesel);   
if (strcmp(cfg.keeptrials,'yes'))
  psth.trial  = singleTrials;
  psth.dimord = 'rpt_chan_time'; % perhaps we need two dimords here, what's the FT convention?
else
  psth.dimord = 'chan_time';
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
if isfield(spike,'cfg'),cfg.previous = spike.cfg; end
% remember the exact configuration details in the output 
psth.cfg     = cfg;


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
  error('ft:spike_psth:cfg:latency:unknownOption',...
  'cfg.latency should be "maxperiod", "minperiod", "prestim", "poststim" or 1-by-2 numerical vector');
end
if cfg.latency(1)>=cfg.latency(2), % check if it is ordered vector
   error('ft:spike_psth:cfg:latency:decreasingLatencyVector',...
   'cfg.latency should be ascending vector, such that cfg.latency(2)>cfg.latency(1)');
end
% check whether the time window fits with the data
if (cfg.latency(1) < min(begTrialLatency)), cfg.latency(1) = min(begTrialLatency); 
  warning('ft:spike_psth:cfg:latencyBegBeforeTrialBeg',...
  'Correcting begin latency because it is before all trial beginnings');
end
if (cfg.latency(2) > max(endTrialLatency)), cfg.latency(2) = max(endTrialLatency);
  warning('ft:spike_psth:cfg:latencyEndAfterTrialEnd',...
  'Correcting end latency because it is after all trial ends');
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
    error('ft:spike_psth:cfg:trials:unknownOption',...
    'cfg.trials should be logical or numerical selection or string "all"');  
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>nTrials, error('ft:spike_psth:cfg:trials:maxExceeded',...
   'maximum trial number in cfg.trials should not exceed length of DATA.trial')
end
if isempty(cfg.trials), error('ft:spike_psth:cfg:trials:noneSelected',...
   'No trials were selected by you, rien ne va plus'); 
end



