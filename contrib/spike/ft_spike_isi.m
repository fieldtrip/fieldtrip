function [isih] = ft_spike_isi(cfg, spike)

% FT_SPIKE_ISI computes the interspike interval histogram
%
% The input SPIKE should be organised as 
% a) the spike datatype, obtained from FT_SPIKE_MAKETRIALS 
% b) the raw datatype, containing binary spike trains, obtained from
% FT_APPENDSPIKE or FT_CHECKDATA. In this case the raw datatype is
% converted to the spike datatype.
%
% Use as
%   [isih] = ft_spike_isi(cfg, spike)
%
% Configurations:
%   cfg.outputunit       = 'spikecount' (default) or 'proportion' (sum of all bins = 1).
%   cfg.spikechannel     = string or index of spike channels to
%                          trigger on (default = 'all')
%                          See FT_CHANNELSELECTION for details.
%   cfg.trials           = numeric selection of trials (default = 'all')
%   cfg.bins             = ascending vector of isi bin edges.
%   cfg.latency          = [begin end] in seconds, 'max' (default), 'min', 'prestim'(t<=0), or
%                          'poststim' (t>=0).
%                          If 'max', we use all available latencies.
%                          If 'min', we use only the time window contained by all trials.
%                          If 'prestim' or 'poststim', we use time to or
%                          from 0, respectively.
%`  cfg.keeptrials       = 'yes' or 'no'. If 'yes', we keep the individual
%                           isis between spikes and output as isih.isi
%   cfg.param            = string, one of
%      'gamfit'      : returns [shape scale] for gamma distribution fit
%      'coeffvar'    : coefficient of variation (sd / mean)
%      'lv'          : Shinomoto's Local Variation measure (2009)
%
% Outputs:
%   isih.avg             = nUnits-by-nBins interspike interval histogram
%   isih.time            = 1 x nBins bincenters corresponding to isih.avg
%   isih.isi             = 1-by-nUnits cell with interval to previous spike per spike.
%                          For example isih.isi{1}(2) = 0.1 means that the
%                          second spike fired was 0.1 s later than the
%                          first. Note that jumps within trials or first
%                          spikes within trials are given NaNs.
%   isih.label           = 1-by-nUnits cell-array with labels

% Copyright (C) 2010-2012, Martin Vinck
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

% check if data is of proper format 
spike = ft_checkdata(spike,'datatype', 'spike', 'feedback', 'yes');

% get the default options
cfg.outputunit   = ft_getopt(cfg,'outputunit','spikecount');
cfg.spikechannel = ft_getopt(cfg,'spikechannel', 'all');
cfg.trials       = ft_getopt(cfg,'trials', 'all');
cfg.latency      = ft_getopt(cfg,'latency','maxperiod');
cfg.keeptrials   = ft_getopt(cfg,'keeptrials', 'yes');
cfg.bins         = ft_getopt(cfg,'bins', 0:0.002:1);
cfg.param        = ft_getopt(cfg,'param', 'coeffvar');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'outputunit','char', {'spikecount', 'proportion'});
cfg = ft_checkopt(cfg,'bins', 'ascendingdoublevector');
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'latency', {'char', 'ascendingdoublebivector'});
cfg = ft_checkopt(cfg,'trials', {'char', 'doublevector', 'logical'}); 
cfg = ft_checkopt(cfg,'keeptrials', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'param', 'char', {'gamfit', 'coeffvar', 'lv'});

cfg = ft_checkconfig(cfg, 'allowed', {'param', 'outputunit', 'bins', 'spikechannel', 'latency', 'trials', 'keeptrials'});

% get the number of trials or change DATA according to cfg.trials
if  strcmp(cfg.trials,'all')
  cfg.trials = 1:size(spike.trialtime,1);
elseif islogical(cfg.trials) || all(cfg.trials==0 | cfg.trials==1)
  cfg.trials = find(cfg.trials);
end
cfg.trials = sort(cfg.trials);

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

% construct the isi bins
bins  = cfg.bins;
nBins = length(bins)-1;

% compute the interspike interval
keepTrials = strcmp(cfg.keeptrials,'yes');
if keepTrials, isiSpike = cell(1,nUnits); end  % contains the individual histogram spike times
isihist  = zeros(nUnits,nBins+1); % isi histogram
out = [];
for iUnit = 1:nUnits
  unitIndx = spikesel(iUnit);
  
  % only select the spikes in the right latencies and trials
  spikeTrials    = spike.trial{unitIndx}(:)'; % ensure row vector
  spikeTimes     = spike.time{unitIndx}(:)'; % ensure row vector
  
  % select only the spikes in the window and with the selected trials
  spikesInTrials = ismember(spikeTrials, cfg.trials); % row vec
  spikesInWin    = spikeTimes>=cfg.latency(1)&spikeTimes<=cfg.latency(2); % row vec
  spikeTimes     = spikeTimes(spikesInTrials & spikesInWin);
  spikeTrials    = spikeTrials(spikesInTrials & spikesInWin);
   
  % find the spikes that jumped to the next trial, replace them with NaNs
  trialJump = logical([1 diff(spikeTrials(:)')]);
  isi = [NaN diff(spikeTimes(:)')];
  isi(trialJump) = NaN;
  
  switch cfg.param
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

  % convert and store the isi
  if keepTrials, isiSpike{iUnit} = isi; end
  isihist(iUnit,:)   = histc(isi,bins);
end

isihist(:,end) = []; % the last number is only an equality to a bin edge
if strcmp(cfg.outputunit,'proportion'),
  isihist = isihist./repmat(nansum(isihist,2),1,size(isihist,2));
end

% gather the rest of the results
if strcmp(cfg.keeptrials,'yes'), isih.isi  = isiSpike; end
isih.time        = bins(1:end-1);
isih.avg         = isihist;
isih.dimord      = 'chan_time';
isih.label       = spike.label(spikesel);
param = cfg.param;
isih.(param) = out;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble previous   spike
ft_postamble provenance isih
ft_postamble history    isih


