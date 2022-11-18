function [psth] = ft_spike_psth(cfg, spike)

% FT_SPIKE_PSTH computes the peristimulus histogram of spiketrains.
%
% Use as
%   [psth] = ft_spike_psth(cfg, spike)
%
% The input SPIKE should be organised as either the spike datatype,
% obtained from FT_SPIKE_MAKETRIALS, or the raw datatype, containing binary
% spike trains, obtained from FT_APPENDSPIKE or FT_CHECKDATA. In this case
% the raw datatype is converted to the spike datatype.
%
% Configurations:
%   cfg.binsize          =  [binsize] in sec or string.
%                          If 'scott', we estimate the optimal bin width
%                          using Scott's formula (1979). If 'sqrt', we take
%                          the number of bins as the square root of the
%                          number of observations. The optimal bin width is
%                          derived over all neurons; thus, this procedure
%                          works best if the input contains only one neuron
%                          at a time.
%   cfg.outputunit       = 'rate' (default) or 'spikecount' or
%                          'proportion'. If 'rate', we
%                          convert the output per trial to firing rates
%                          (spikes/sec). If 'spikecount', we count the
%                          number spikes per trial. If 'proportion', we
%                          normalize the area under the PSTH to 1.
%   cfg.spikechannel     = See FT_CHANNELSELECTION for details. cfg.trials
%                          is vector of indices (e.g., 1:2:10)
%                          logical selection of trials (e.g., [1010101010])
%                          'all' (default), selects all trials
%   cfg.vartriallen      = 'yes' (default)
%                          Accept variable trial lengths and use all
%                          available trials and the samples in every trial.
%                          Missing values will be ignored in the
%                          computation of the average and the variance and
%                          stored as NaNs in the output psth.trial. 'no'
%                          Only select those trials that fully cover the
%                          window as specified by cfg.latency and discard
%                          those trials that do not.
%   cfg.latency          = [begin end] in seconds
%                          'maxperiod' (default), i.e., maximum period
%                          available 'minperiod', i.e., the minimal period
%                          all trials share, 'prestim' (all t<=0) 'poststim'
%                          (all t>=0).
%   cfg.keeptrials       = 'yes' or 'no' (default)
%   cfg.trials           =  numeric or logical selection of trials (default = 'all')
%
% Outputs:
%   Psth is a timelock datatype (see FT_DATATYPE_TIMELOCK)
%     Psth.time        = center histogram bin points
%	    Psth.fsample     = 1/binsize;
%     Psth.avg         = contains average PSTH per unit
%     Psth.trial       = contains PSTH per unit per trial
%     Psth.var         = contains variance of PSTH per unit across trials
%
% For subsequent processing you can use
%   FT_SPIKE_PLOT_PSTH    : plot only the PSTH, for a single neuron
%   FT_TIMELOCKSTATISTICS : compute statistics on the PSTH
%   FT_SPIKE_PLOT_RASTER  : plot PSTH with raster for one or more neurons
%   FT_SPIKE_JPSTH        : compute the JPSTH

%  Copyright (C) 2010-2013, Martin Vinck
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


% control input spike structure and convert to spike if raw structure
spike = ft_checkdata(spike,'datatype', 'spike', 'feedback', 'yes');

% get the default options
cfg.outputunit   = ft_getopt(cfg, 'outputunit','rate');
cfg.binsize      = ft_getopt(cfg, 'binsize','scott');
cfg.spikechannel = ft_getopt(cfg, 'spikechannel', 'all');
cfg.trials       = ft_getopt(cfg, 'trials', 'all');
cfg.latency      = ft_getopt(cfg, 'latency','maxperiod');
cfg.vartriallen  = ft_getopt(cfg, 'vartriallen', 'yes');
cfg.keeptrials   = ft_getopt(cfg, 'keeptrials', 'no');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'outputunit',   'char',  {'rate', 'spikecount'});
cfg = ft_checkopt(cfg,'binsize',     {'char',  'doublescalar'});
cfg = ft_checkopt(cfg,'spikechannel',{'cell',  'char', 'double'});
cfg = ft_checkopt(cfg,'latency',     {'char',  'ascendingdoublebivector'});
cfg = ft_checkopt(cfg,'trials',      {'char',  'doublevector', 'logical'});
cfg = ft_checkopt(cfg,'vartriallen' , 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'keeptrials'  , 'char', {'yes', 'no'});

cfg = ft_checkconfig(cfg, 'allowed', {'outputunit', 'binsize', 'spikechannel', 'trials', 'latency', 'vartriallen', 'keeptrials'});

% get the number of trials or convert to indices
cfg        = trialselection(cfg,spike);

% select the unit - this should be done with channelselection function
cfg.spikechannel = ft_channelselection(cfg.spikechannel, spike.label);
spikesel         = match_str(spike.label, cfg.spikechannel);
nUnits           = length(spikesel);
if nUnits==0, error('no spikechannel selected by means of cfg.spikechannel'); end

% determine the duration of each trial
begTrialLatency = spike.trialtime(cfg.trials,1); % remember: already selected on trial here
endTrialLatency = spike.trialtime(cfg.trials,2);
trialDur 	    	= endTrialLatency - begTrialLatency;

% select the latencies, use the same modular function in all the scripts
cfg = latencyselection(cfg,begTrialLatency,endTrialLatency);

% compute the optimal bin width if desired
if ischar(cfg.binsize)
  h = zeros(1,nUnits);
  for iUnit = 1:nUnits
    unitIndx       = spikesel(iUnit); % select the unit
    spikesInWin    = spike.time{unitIndx}>=cfg.latency(1) & spike.time{unitIndx}<=cfg.latency(2); % get spikes in trial
    
    % automatically determine an 'optimal' binwidth
    N              = sum(spikesInWin);
    if strcmp(cfg.binsize,'scott')
      sd           = nanstd(spike.time{unitIndx}(spikesInWin));
      h(iUnit)     = 3.49*sd./(N^(1/3));
    elseif strcmp(cfg.binsize,'sqrt')
      k            = ceil(sqrt(N));
      h(iUnit)     = (cfg.latency(2)-cfg.latency(1))/k;
    else
      error('unsupported option for cfg.binsize');
    end
  end
  cfg.binsize = nanmean(h);
end

% do some error checking on the binsize
if cfg.binsize<=0 || cfg.binsize>(cfg.latency(2)-cfg.latency(1))
  error('cfg.binsize should be greater than zero and not exceed the trialduration');
end

% end of trial should be late enough, beginning should be early enough
fullDur       = trialDur>=cfg.binsize; % trials which have full duration
overlaps      = endTrialLatency>=(cfg.latency(1)+cfg.binsize) & begTrialLatency<=(cfg.latency(2)-cfg.binsize);
if strcmp(cfg.vartriallen,'no') % only select trials that fully cover our latency window
  startsLater    = single(begTrialLatency) > (single(cfg.latency(1)) + 0.0001); % last factor just for rounding errors
  endsEarlier    = single(endTrialLatency) < (single(cfg.latency(2)) - 0.0001);
  hasWindow      = ~(startsLater | endsEarlier); % it should not start later or end earlier
else
  hasWindow     = ones(1,length(cfg.trials));
end
trialSel        = fullDur(:) & overlaps(:) & hasWindow(:);
cfg.trials      = cfg.trials(trialSel); % note that endTrialLatency was of length cfg.trials
if isempty(cfg.trials), warning('No trials were selected after latency selection'); end
nTrials         = length(cfg.trials);
begTrialLatency = begTrialLatency(trialSel); % note that begTrialLatency was of length cfg.trials here
endTrialLatency = endTrialLatency(trialSel);

% create the bins, we start arbitrarily at cfg.latency(1)
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

for iTrial = 1:nTrials
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
    if strcmp(cfg.outputunit,'rate'),
      trialPsth = trialPsth/cfg.binsize;
    elseif strcmp(cfg.outputunit,'proportion'),
      trialPsth = trialPsth./nansum(trialPsth);
    end
    
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
psth.time      = bins(1:end-1) + 0.5*cfg.binsize;
psth.label     = spike.label(spikesel);
if (strcmp(cfg.keeptrials,'yes'))
  psth.trial  = singleTrials;
  psth.dimord = 'rpt_chan_time';
else
  psth.dimord = 'chan_time';
end
if isfield(spike,'sampleinfo'), psth.sampleinfo = spike.sampleinfo(cfg.trials,:); end
if isfield(spike,'trialinfo'),  psth.trialinfo  = spike.trialinfo(cfg.trials,:);  end

% do the general cleanup and bookkeeping at the end of the function

ft_postamble previous   spike
ft_postamble provenance psth
ft_postamble history    psth


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
% check whether the time window fits with the data
if (cfg.latency(1) < min(begTrialLatency)), cfg.latency(1) = min(begTrialLatency);
  warning('Correcting begin latency because it is before all trial beginnings');
end
if (cfg.latency(2) > max(endTrialLatency)), cfg.latency(2) = max(endTrialLatency);
  warning('Correcting end latency because it is after all trial ends');
end


function [cfg] = trialselection(cfg,spike)

nTrials = size(spike.trialtime,1);
if  strcmp(cfg.trials,'all')
  cfg.trials = 1:nTrials;
elseif islogical(cfg.trials) || all(cfg.trials==0 | cfg.trials==1)
  cfg.trials = find(cfg.trials);
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>nTrials,
  error('maximum trial number in cfg.trials should not exceed number of rows of spike.trialtime');
end
if isempty(cfg.trials), error('No trials were selected by you, rien ne va plus'); end
