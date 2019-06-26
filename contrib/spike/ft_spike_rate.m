function [rate] = ft_spike_rate(cfg, spike)

% FT_SPIKE_RATE computes the firing rate of spiketrains and their variance
%
% Use as
%   [rate] = ft_spike_rate(cfg, spike)
%
% The input SPIKE should be organised as the spike or the (binary) raw
% datatype, obtained from FT_SPIKE_MAKETRIALS or FT_APPENDSPIKE (in that
% case, conversion is done within the function)
%
% Configurations:
%   cfg.outputunit       = 'rate' (default) or 'spikecount'. If 'rate', we convert
%                          the output per trial to firing rates (spikes/sec).
%                          If 'spikecount', we count the number spikes per trial.
%   cfg.spikechannel     = see FT_CHANNELSELECTION for details
%   cfg.trials           = vector of indices (e.g., 1:2:10)
%                          logical selection of trials (e.g., [1010101010])
%                          'all' (default), selects all trials%   cfg.trials
%   cfg.vartriallen      = 'yes' (default) or 'no'.
%                          If 'yes' - accept variable trial lengths and use all available trials
%                          and the samples in every trial.
%                          If 'no'  - only select those trials that fully cover the window as
%                          specified by cfg.latency and discard those trials that do not.
%   cfg.latency          = [begin end] in seconds
%                          'maxperiod' (default)
%                          'minperiod', i.e., the minimal period all trials share
%                          'prestim' (all t<=0)
%                          'poststim' (all t>=0).
%   cfg.keeptrials       = 'yes' or 'no' (default).
%
% The outputs from spike is a TIMELOCK structure (see FT_DATATYPE_TIMELOCK)
%   rate.trial:        nTrials x nUnits matrix containing the firing rate per unit and trial
%   rate.avg:          nTrials array containing the average firing rate per unit
%   rate.var:          nTrials array containing the variance of firing rates per unit
%   rate.dof:          nTrials array containing the degree of freedom per unit
%   rate.label:        nUnits cell-array containing the labels of the neuronal units%
%   rate.time:         Mean latency (this field ensures it is TIMELOCK
%                      struct)

% Copyright (C) 2010, Martin Vinck
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
ft_preamble provenance spike
ft_preamble trackconfig

% control input spike structure
spike = ft_checkdata(spike,'datatype', 'spike', 'feedback', 'yes');

% get the default options
if isfield(cfg,'trials') && isempty(cfg.trials), error('no trials were selected'); end % empty should result in error, not in default
cfg.outputunit   = ft_getopt(cfg, 'outputunit','rate');
cfg.spikechannel = ft_getopt(cfg, 'spikechannel', 'all');
cfg.trials       = ft_getopt(cfg, 'trials', 'all');
cfg.latency      = ft_getopt(cfg,'latency','maxperiod');
cfg.vartriallen  = ft_getopt(cfg,'vartriallen', 'yes');
cfg.keeptrials   = ft_getopt(cfg,'keeptrials', 'yes');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'outputunit','char', {'rate', 'spikecount'});
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'latency', {'char', 'ascendingdoublebivector'});
cfg = ft_checkopt(cfg,'trials', {'char', 'doublevector', 'logical'}); 
cfg = ft_checkopt(cfg,'vartriallen', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'keeptrials', 'char', {'yes', 'no'});

% check if the configuration inputs are valid
cfg = ft_checkconfig(cfg, 'allowed', {'outputunit', 'spikechannel', 'trials', 'latency', 'vartriallen', 'keeptrials'});

% get the spikechannels
cfg.spikechannel = ft_channelselection(cfg.spikechannel, spike.label);
spikesel    = match_str(spike.label, cfg.spikechannel);
nUnits      = length(spikesel); % number of spike channels
if nUnits==0, error('No spikechannel selected by means of cfg.spikechannel'); end

% do the trial selection
cfg        = trialselection(cfg,spike);
nTrials    = length(cfg.trials); % actual number of trials we use

% determine the duration of each trial
begTrialLatency = spike.trialtime(cfg.trials,1);
endTrialLatency = spike.trialtime(cfg.trials,2);

% select the latencies
cfg = latencyselection(cfg,begTrialLatency,endTrialLatency);

% check which trials will be used based on the latency
overlaps      = endTrialLatency>(cfg.latency(1)) & begTrialLatency<(cfg.latency(2));
hasWindow     = true(nTrials,1);
if strcmp(cfg.vartriallen,'no') % only select trials that fully cover our latency window
  startsLater    = begTrialLatency > (cfg.latency(1) + 0.0001);
  endsEarlier    = endTrialLatency < (cfg.latency(2) - 0.0001);
  hasWindow      = ~(startsLater | endsEarlier); % it should not start later or end earlier
  trialDur       = ones(nTrials,1)*(cfg.latency(2)-cfg.latency(1));
elseif strcmp(cfg.vartriallen,'yes')
  winBeg      = max([begTrialLatency(:) cfg.latency(1)*ones(nTrials,1)],[],2); 
  winEnd      = min([endTrialLatency(:) cfg.latency(2)*ones(nTrials,1)],[],2);
  trialDur    = winEnd-winBeg; % the effective trial duration  
end
cfg.trials         = cfg.trials(overlaps(:) & hasWindow(:));
trialDur           = trialDur(overlaps(:) & hasWindow(:)); % select the durations, we will need this later
nTrials            = length(cfg.trials);
if isempty(cfg.trials), warning('No trials were selected in the end'); end

% preallocate before computing the psth
keepTrials = strcmp(cfg.keeptrials,'yes');
if  keepTrials, singleTrials = NaN(nTrials,nUnits); end % preallocate single trials with NaNs
[s,ss]   = deal(zeros(nUnits,1));
dof      = nTrials*ones(nUnits,1); % compute degrees of freedom
for iUnit = 1:nUnits
  unitIndx   = spikesel(iUnit);
  ts         = spike.time{unitIndx}(:); % get the times
  latencySel = ts>=cfg.latency(1) & ts<=cfg.latency(2); % select timestamps within latency
  trialSel   = ismember(spike.trial{unitIndx},cfg.trials);
  trialNums  = spike.trial{unitIndx}(latencySel(:) & trialSel(:));
  
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
rate.dimord    = 'chan_time';
rate.time      = mean(cfg.latency);
if (strcmp(cfg.keeptrials,'yes'))
  rate.trial = singleTrials;
  rate.dimord = 'rpt_chan_time';
end
if isfield(spike, 'trialinfo'), rate.trialinfo = spike.trialinfo(cfg.trials,:); end
  
% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble previous   spike
ft_postamble provenance rate
ft_postamble history    rate


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
end

function [cfg] = trialselection(cfg,spike)

% get the number of trials or change DATA according to cfg.trials
nTrials = size(spike.trialtime,1);
if  strcmp(cfg.trials,'all')
  cfg.trials = 1:nTrials;
elseif islogical(cfg.trials) || all(cfg.trials==0 | cfg.trials==1)
  cfg.trials = find(cfg.trials);
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>nTrials, error('maximum trial number in cfg.trials should not exceed size(spike.trialtime,1)'); end
if isempty(cfg.trials), error('no trials were selected by you'); end

