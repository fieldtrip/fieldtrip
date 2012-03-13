function [timelock] = ft_spiketriggeredaverage(cfg, data)

% FT_SPIKETRIGGEREDAVERAGE computes the avererage of the LFP around the
% spikes.
%
% Use as
%   [timelock] = ft_spiketriggeredaverage(cfg, data)
%
% The input data should be organised in a structure as obtained from
% the FT_PREPROCESSING function. The configuration should be according to
%
%   cfg.timwin       = [begin end], time around each spike (default = [-0.1 0.1])
%   cfg.spikechannel = string, name of single spike channel to trigger on
%   cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                      see FT_CHANNELSELECTION for details
%   cfg.keeptrials   = 'yes' or 'no', return individual trials or average (default = 'no')
%   cfg.feedback     = 'no', 'text', 'textbar', 'gui' (default = 'no')

% Copyright (C) 2008, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% check input data structure
data = ft_checkdata(data,'datatype', 'raw', 'feedback', 'yes');

% these were supported in the past, but are not any more (for consistency with other spike functions)
cfg = ft_checkconfig(cfg, 'forbidden', 'inputfile', ...
                                       'outputfile');  % see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1056

%get the options
cfg.timwin       = ft_getopt(cfg, 'timwin',[-0.1 0.1]);
cfg.spikechannel = ft_getopt(cfg,'spikechannel', []);
cfg.channel      = ft_getopt(cfg,'channel', 'all');
cfg.keeptrials   = ft_checkopt(cfg,'keeptrials', 'char', {'yes', 'no'});
cfg.feedback     = ft_checkopt(cfg,'feedback', 'yes');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'timwin','doublevector');
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'channel', {'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'keeptrials', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'feedback', 'char', {'yes', 'no'});

% autodetect the spike channels
ntrial = length(data.trial);
nchans  = length(data.label);
spikechan = zeros(nchans,1);
for i=1:ntrial
  for j=1:nchans
    spikechan(j) = spikechan(j) + all(data.trial{i}(j,:)==0 | data.trial{i}(j,:)==1 | data.trial{i}(j,:)==2);
  end
end
spikechan = (spikechan==ntrial);

% determine the channels to be averaged
cfg.channel = ft_channelselection(cfg.channel, data.label);
chansel     = match_str(data.label, cfg.channel);
nchansel    = length(cfg.channel);  % number of channels

% determine the spike channel on which will be triggered
cfg.spikechannel = ft_channelselection(cfg.spikechannel, data.label);
spikesel         = match_str(data.label, cfg.spikechannel);
nspikesel        = length(cfg.spikechannel);    % number of channels

if nspikesel==0
  error('no spike channel selected');
end

if nspikesel>1
  error('only supported for a single spike channel');
end

if ~spikechan(spikesel)
  error('the selected spike channel seems to contain continuous data');
end

begpad = round(cfg.timwin(1)*data.fsample);
endpad = round(cfg.timwin(2)*data.fsample);
numsmp = endpad - begpad + 1;

singletrial = cell(1,ntrial);
spiketime   = cell(1,ntrial);
spiketrial  = cell(1,ntrial);
cumsum = zeros(nchansel, numsmp);
cumcnt = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:ntrial
  spikesmp = find(data.trial{i}(spikesel,:));
  spikecnt = data.trial{i}(spikesel,spikesmp);
  
  if any(spikecnt>5) || any(spikecnt<0)
    error('the spike count lies out of the regular bounds');
  end
  
  % instead of doing the bookkeeping of double spikes below, replicate the double spikes by looking at spikecnt
  sel = find(spikecnt>1);
  tmp = zeros(1,sum(spikecnt(sel)));
  n   = 1;
  for j=1:length(sel)
    for k=1:spikecnt(sel(j))
      tmp(n) = spikesmp(sel(j));
      n = n + 1;
    end
  end
  spikesmp(sel) = [];                     % remove the double spikes
  spikecnt(sel) = [];                     % remove the double spikes
  spikesmp = [spikesmp tmp];              % add the double spikes as replicated single spikes
  spikecnt = [spikecnt ones(size(tmp))];  % add the double spikes as replicated single spikes
  spikesmp = sort(spikesmp);              % sort them to keep the original ordering (not needed on spikecnt, since that is all ones)
  
  spiketime{i}  = data.time{i}(spikesmp);
  spiketrial{i} = i*ones(size(spikesmp));
  fprintf('processing trial %d of %d (%d spikes)\n', i, ntrial, sum(spikecnt));
  
  if strcmp(cfg.keeptrials, 'yes')
    if any(spikecnt>1)
      error('overlapping spikes not supported with cfg.keeptrials=yes');
    end
    % initialize the memory for this trial
    singletrial{i} = nan*zeros(length(spikesmp), nchansel, numsmp);
  end
  
  ft_progress('init', cfg.feedback, 'averaging spikes');
  for j=1:length(spikesmp)
    ft_progress(i/ntrial, 'averaging spike %d of %d\n', j, length(spikesmp));
    begsmp = spikesmp(j) + begpad;
    endsmp = spikesmp(j) + endpad;
    
    if begsmp<1
      % a possible alternative would be to pad the begin with nan
      % this excludes the complete segment
      continue
    elseif endsmp>size(data.trial{i},2)
      % possible alternative would be to pad the end with nan
      % this excludes the complete segment
      continue
    else
      segment = data.trial{i}(chansel,begsmp:endsmp);
    end
    if strcmp(cfg.keeptrials, 'yes')
      singletrial{i}(j,:,:) = segment;
    end
    
    cumsum = cumsum + spikecnt(j)*segment;
    cumcnt = cumcnt + spikecnt(j);
    
  end % for each spike in this trial
  ft_progress('close');
  
end % for each trial


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timelock.time  = offset2time(begpad, data.fsample, numsmp);
timelock.avg   = cumsum ./ cumcnt;
timelock.label = data.label(chansel);

if (strcmp(cfg.keeptrials, 'yes'))
  timelock.dimord = 'rpt_chan_time';
  % concatenate all the single spike snippets
  timelock.trial     = cat(1, singletrial{:});
  timelock.origtime  = cat(2,spiketime{:})';  % this deviates from the standard output, but is included for reference
  timelock.origtrial = cat(2,spiketrial{:})'; % this deviates from the standard output, but is included for reference
  
  % select all trials that do not contain data in the first sample
  sel = isnan(timelock.trial(:,1,1));
  fprintf('removing %d trials from the output that do not contain data\n', sum(sel));
  % remove the selected trials from the output
  timelock.trial       = timelock.trial(~sel,:,:);
  timelock.origtime    = timelock.origtime(~sel);
  timelock.origtrial   = timelock.origtrial(~sel);
else
  timelock.dimord = 'chan_time';
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous data
ft_postamble history timelock
