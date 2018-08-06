function [spike] = ft_spike_select(cfg, spike)

% FT_SPIKE_SELECT selects subsets of spikes, channels and trials from a
% spike structure.
%
% Use as
%   [spike] = ft_spike_select(cfg, spike)
%
% The input SPIKE should be organised as the spike datatype (see
% FT_DATATYPE_SPIKE) 
%
% Configurations:
%   cfg.spikechannel     = See FT_CHANNELSELECTION for details.
%   cfg.trials           = vector of indices (e.g., 1:2:10)
%                          logical selection of trials (e.g., [1010101010])
%                          'all' (default), selects all trials
%   cfg.latency          = [begin end] in seconds
%                          'maxperiod' (default), i.e., maximum period available
%                          'minperiod', i.e., the minimal period all trials share
%                          'prestim' (all t<=0)
%                          'poststim' (all t>=0).
% Outputs:
%   Spike structure with selections

% Copyright (C) 2012, Martin Vinck
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
cfg.spikechannel = ft_getopt(cfg, 'spikechannel', 'all');
cfg.trials       = ft_getopt(cfg, 'trials', 'all');
cfg.latency      = ft_getopt(cfg, 'latency', 'maxperiod');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'latency', {'char', 'ascendingdoublebivector'});
cfg = ft_checkopt(cfg,'trials', {'char', 'doublevector', 'logical'}); 

% select the desired channels
doAll = strcmp(cfg.spikechannel,'all');
try
  cfg.spikechannel = ft_channelselection(cfg.spikechannel, spike.label);
catch
  [I,J] = unique(spike.label, 'first');
  label = spike.label(J);
  cfg.spikechannel = ft_channelselection(cfg.spikechannel, label);
end  
spikesel         = match_str(spike.label, cfg.spikechannel);
nUnits           = length(spikesel);
if nUnits==0, error('no spikechannel selected by means of cfg.spikechannel');end
doAllTrials = strcmp(cfg.trials,'all'); 
doAllLatencies = strcmp(cfg.latency,'maxperiod'); 

% select the desired channels
if ~doAll
  fprintf('Selecting channels\n');
  try, spike.time = spike.time(spikesel); end
  try, spike.trial = spike.trial(spikesel); end
  try, spike.waveform = spike.waveform(spikesel); end
  try, spike.timestamp = spike.timestamp(spikesel); end
  try, spike.unit      = spike.unit(spikesel); end
  try, spike.fourierspctrm = spike.fourierspctrm(spikesel); end
  try, spike.label = spike.label(spikesel); end
end

% select the desired trials
if ~isfield(spike,'trial') || ~isfield(spike,'trialtime') || ~isfield(spike,'time')
  if ~doAllTrials
    warning('spike structure does not contain trial, time or trialtime field, cannot select trials');
  end
else
  doSelection = ~strcmp(cfg.trials,'all');
  cfg  = trialselection(cfg,spike);  
  if doSelection
    fprintf('Selecting trials\n');      
    nUnits  = length(spike.trial);
    for iUnit = 1:nUnits
      spikesInTrial = ismember(spike.trial{iUnit},cfg.trials); % get spikes in trial
      spike.trial{iUnit} = spike.trial{iUnit}(spikesInTrial);
      spike.time{iUnit} = spike.time{iUnit}(spikesInTrial);
      try, spike.timestamp{iUnit} = spike.timestamp{iUnit}(spikesInTrial); end
      try, spike.waveform{iUnit} = spike.waveform{iUnit}(spikesInTrial); end
      try, spike.unit{iUnit} = spike.unit{iUnit}(spikesInTrial); end
      try, spike.fourierspctrm{iUnit} = spike.fourierspctrm{iUnit}(spikesInTrial,:,:); end
    end
    warning('No selection is performed on spike.trialtime and spike.sampleinfo'); 
    % Note: .trialtime and .sampleinfo are explicitly left intact, as the
    % relationship between trial number and trialtime becomes ambigious
    % otherwise!
  end
end

% select the desired latencies
if ~isfield(spike, 'trialtime') || ~isfield(spike,'time') || ~isfield(spike,'trial')
  if ~doAllLatencies
    warning('cannot select latencies as either .trialtime, .time or .trial is missing'); 
  end
else
  if ~strcmp(cfg.latency,'maxperiod') % otherwise, all spikes are selected
    fprintf('Selecting on latencies\n');
    begTrialLatency = spike.trialtime(cfg.trials,1); 
    endTrialLatency = spike.trialtime(cfg.trials,2);
    cfg = latencyselection(cfg,begTrialLatency,endTrialLatency);  
    for iUnit = 1:nUnits
        spikesInWin = spike.time{iUnit}>=cfg.latency(1) & spike.time{iUnit}<=cfg.latency(2); % get spikes in trial
        spike.trial{iUnit} = spike.trial{iUnit}(spikesInWin);
        spike.time{iUnit} = spike.time{iUnit}(spikesInWin);
        try, spike.timestamp{iUnit} = spike.timestamp{iUnit}(spikesInWin); end
        try, spike.waveform{iUnit} = spike.waveform{iUnit}(spikesInWin); end
        try, spike.unit{iUnit} = spike.unit{iUnit}(spikesInWin); end
        try, spike.fourierspctrm{iUnit} = spike.fourierspctrm{iUnit}(spikesInWin,:,:); end
    end
    % modify the trialtime
    adjust = spike.trialtime(:,1)<=cfg.latency(1);
    spike.trialtime(adjust,1) = cfg.latency(1);
    adjust = spike.trialtime(:,2)>=cfg.latency(2);
    spike.trialtime(adjust,2) = cfg.latency(2);   
    
    % FIXME: change sampleinfo field automatically again, for now remove.    
    try, spike = rmfield(spike,'sampleinfo'); end
  end
end
    
% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble previous   spike
ft_postamble provenance spike
ft_postamble history    spike


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

%%%
function [cfg] = trialselection(cfg,spike)

% get the number of trials or change DATA according to cfg.trials
nTrials = size(spike.trialtime,1);
if  strcmp(cfg.trials,'all')
  cfg.trials = 1:nTrials;
elseif islogical(cfg.trials)
  cfg.trials = find(cfg.trials);
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>nTrials, error('maximum trial number in cfg.trials should not exceed length of spike.trialtime')
end
if isempty(cfg.trials), error('No trials were selected');
end



