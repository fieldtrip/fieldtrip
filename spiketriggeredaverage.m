function [timelock] = spiketriggeredaverage(cfg, data)

% SPIKETRIGGEREDAVERAGE computes the avererage of teh LFP around the spikes.
%
% Use as
%   [timelock] = spiketriggeredaverage(cfg, data)
%
% The input data should be organised in a structure as obtained from
% the PREPROCESSING function. The configuration should be according to
%
%   cfg.timwin       = [begin end], time around each spike (default = [-0.1 0.1])
%   cfg.spikechannel = string, name of single spike channel to trigger on
%   cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                      see CHANNELSELECTION for details
%   cfg.keeptrials   = 'yes' or 'no', return individual trials or average (default = 'no')
%   cfg.feedback     = 'no', 'text', 'textbar', 'gui' (default = 'no')

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: spiketriggeredaverage.m,v $
% Revision 1.5  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.4  2008/04/09 14:39:22  roboos
% remove trials without data (too close to trial border) from the output
%
% Revision 1.3  2008/04/09 14:20:31  roboos
% replicate double spikes, for keeptrials
% give error if spike count >5
%
% Revision 1.2  2008/03/18 21:56:56  roboos
% fixed bug in keeptrials, nans were allocated at the wrong moment
%
% Revision 1.1  2008/03/17 15:08:16  roboos
% new implementation, based on discussion with Thilo about desired functionality replication of spikeanalysis
%

fieldtripdefs

% set the defaults
if ~isfield(cfg, 'timwin'),       cfg.timwin = [-0.1 0.1];    end
if ~isfield(cfg, 'channel'),      cfg.channel = 'all';        end
if ~isfield(cfg, 'spikechannel'), cfg.spikechannel = [];      end
if ~isfield(cfg, 'keeptrials'),   cfg.keeptrials = 'no';      end
if ~isfield(cfg, 'feedback'),     cfg.feedback = 'no';        end

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
cfg.channel = channelselection(cfg.channel, data.label);
chansel     = match_str(data.label, cfg.channel);
nchansel    = length(cfg.channel);	% number of channels

% determine the spike channel on which will be triggered
cfg.spikechannel = channelselection(cfg.spikechannel, data.label);
spikesel         = match_str(data.label, cfg.spikechannel);
nspikesel        = length(cfg.spikechannel);	% number of channels

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

  progress('init', cfg.feedback, 'averaging spikes');
  for j=1:length(spikesmp)
    progress(i/ntrial, 'averaging spike %d of %d\n', j, length(spikesmp));
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
  progress('close');

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

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: spiketriggeredaverage.m,v 1.5 2008/09/22 20:17:44 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
timelock.cfg = cfg;


