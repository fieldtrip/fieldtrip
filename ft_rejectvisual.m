function [data] = ft_rejectvisual(cfg, data)

% FT_REJECTVISUAL shows the preprocessed data in all channels and/or trials to
% allow the user to make a visual selection of the data that should be
% rejected. The data can be displayed in a "summary" mode, in which case
% the variance (or another metric) in each channel and each trial is
% computed. Alternatively, all channels can be shown at once allowing
% paging through the trials, or all trials can be shown, allowing paging
% through the channels.
%
% Use as
%   [data] = ft_rejectvisual(cfg, data)
%
% The configuration can contain
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'),
%                     see FT_CHANNELSELECTION for details
%   cfg.trials      = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.latency     = [begin end] in seconds, or 'minperlength', 'maxperlength',
%                     'prestim', 'poststim' (default = 'maxperlength')
%   cfg.method      = string, describes how the data should be shown, this can be
%                     'summary'  show a single number for each channel and trial (default)
%                     'channel'  show the data per channel, all trials at once
%                     'trial'    show the data per trial, all channels at once
%   cfg.keepchannel = string, determines how to deal with channels that are
%                     not selected, can be
%                     'no'     completely remove unselected channels from the data (default)
%                     'yes'    keep unselected channels in the output data
%                     'nan'    fill the channels that are unselected with NaNs
%   cfg.metric      = string, describes the metric that should be computed in summary mode
%                     for each channel in each trial, can be
%                     'var'       variance within each channel (default)
%                     'min'       minimum value in each channel
%                     'max'       maximum value each channel
%                     'maxabs'    maximum absolute value in each channel
%                     'range'     range from min to max in each channel
%                     'kurtosis'  kurtosis, i.e. measure of peakedness of the amplitude distribution
%                     'zvalue'    mean and std computed over all time and trials, per channel
%   cfg.alim        = value that determines the amplitude scaling for the
%                     channel and trial display, if empty then the amplitude
%                     scaling is automatic (default = [])
%   cfg.eegscale    = number, scaling to apply to the EEG channels prior to display
%   cfg.eogscale    = number, scaling to apply to the EOG channels prior to display
%   cfg.ecgscale    = number, scaling to apply to the ECG channels prior to display
%   cfg.emgscale    = number, scaling to apply to the EMG channels prior to display
%   cfg.megscale    = number, scaling to apply to the MEG channels prior to display
%   cfg.gradscale   = number, scaling to apply to the MEG gradiometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.magscale    = number, scaling to apply to the MEG magnetometer channels prior to display (in addition to the cfg.megscale factor)
%
% The scaling to the EEG, EOG, ECG, EMG and MEG channels is optional and can
% be used to bring the absolute numbers of the different channel types in
% the same range (e.g. fT and uV). The channel types are determined from
% the input data using FT_CHANNELSELECTION.
%
% Optionally, the raw data is preprocessed (filtering etc.) prior to
% displaying it or prior to computing the summary metric. The
% preprocessing and the selection of the latency window is NOT applied
% to the output data.
%
% The following settings are usefull for identifying EOG artifacts:
%   cfg.preproc.bpfilter    = 'yes'
%   cfg.preproc.bpfilttype  = 'but'
%   cfg.preproc.bpfreq      = [1 15]
%   cfg.preproc.bpfiltord   = 4
%   cfg.preproc.rectify     = 'yes'
%
% The following settings are usefull for identifying muscle artifacts:
%   cfg.preproc.bpfilter    = 'yes'
%   cfg.preproc.bpfreq      = [110 140]
%   cfg.preproc.bpfiltord   =  8
%   cfg.preproc.bpfilttype  = 'but'
%   cfg.preproc.rectify     = 'yes'
%   cfg.preproc.boxcar      = 0.2
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_REJECTARTIFACT, FT_REJECTCOMPONENT

% Undocumented local options:
% cfg.feedback

% Copyright (C) 2005-2006, Markus Bauer, Robert Oostenveld
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

% Undocumented options
% cfg.plotlayout = 'square' (default) or '1col', plotting every channel/trial under each other
% cfg.viewmode   = 'remove' (default) or 'toggle', remove the data points from the plot, or mark them (summary mode), which allows for getting them back

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar data

% ft_checkdata is done further down

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval',  {'metric',  'absmax',  'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method',  'absmax',  'maxabs'});

% set the defaults
if ~isfield(cfg, 'channel'),     cfg.channel = 'all';          end
if ~isfield(cfg, 'trials'),      cfg.trials = 'all';           end
if ~isfield(cfg, 'latency'),     cfg.latency = 'maxperlength'; end
if ~isfield(cfg, 'keepchannel'), cfg.keepchannel = 'no';       end
if ~isfield(cfg, 'feedback'),    cfg.feedback = 'textbar';     end
if ~isfield(cfg, 'method'),      cfg.method = 'summary';       end
if ~isfield(cfg, 'alim'),        cfg.alim = [];                end
if ~isfield(cfg, 'eegscale'),    cfg.eegscale = [];            end
if ~isfield(cfg, 'eogscale'),    cfg.eogscale = [];            end
if ~isfield(cfg, 'ecgscale'),    cfg.ecgscale = [];            end
if ~isfield(cfg, 'emgscale'),    cfg.emgscale = [];            end
if ~isfield(cfg, 'megscale'),    cfg.megscale = [];            end
if ~isfield(cfg, 'gradscale'),   cfg.gradscale = [];           end
if ~isfield(cfg, 'magscale'),    cfg.magscale = [];            end
if ~isfield(cfg, 'plotlayout'),  cfg.plotlayout = 'square';    end
if ~isfield(cfg, 'viewmode'),    cfg.viewmode   = 'remove';    end

% store original datatype
dtype = ft_datatype(data);

% check if the input data is valid for this function, this will convert it to raw if needed
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes', 'hassampleinfo', 'yes');

% for backward compatibility
if ~isfield(cfg, 'metric') && any(strcmp(cfg.method, {'var', 'min', 'max', 'maxabs', 'range'}))
  cfg.metric = cfg.method;
  cfg.method = 'summary';
end

if ~isfield(cfg, 'metric')
  cfg.metric = 'var';
end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data = ft_selectdata(data, 'rpt', cfg.trials);
end

% determine the duration of each trial
for i=1:length(data.time)
  begsamplatency(i) = min(data.time{i});
  endsamplatency(i) = max(data.time{i});
end

% determine the latency window which is possible in all trials
minperlength = [max(begsamplatency) min(endsamplatency)];
maxperlength = [min(begsamplatency) max(endsamplatency)];

% latency window for averaging and variance computation is given in seconds
if (strcmp(cfg.latency, 'minperlength'))
  cfg.latency = [];
  cfg.latency(1) = minperlength(1);
  cfg.latency(2) = minperlength(2);
elseif (strcmp(cfg.latency, 'maxperlength'))
  cfg.latency = [];
  cfg.latency(1) = maxperlength(1);
  cfg.latency(2) = maxperlength(2);
elseif (strcmp(cfg.latency, 'prestim'))
  cfg.latency = [];
  cfg.latency(1) = maxperlength(1);
  cfg.latency(2) = 0;
elseif (strcmp(cfg.latency, 'poststim'))
  cfg.latency = [];
  cfg.latency(1) = 0;
  cfg.latency(2) = maxperlength(2);
end

% ensure that the preproc specific options are located in the cfg.preproc substructure
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});

% apply scaling to the selected channel types to equate the absolute numbers (i.e. fT and uV)
% make a seperate copy to prevent the original data from being scaled
tmpdata = data;
scaled  = 0;
if ~isempty(cfg.eegscale)
  scaled = 1;
  chansel = match_str(tmpdata.label, ft_channelselection('EEG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.eegscale;
  end
end
if ~isempty(cfg.eogscale)
  scaled = 1;
  chansel = match_str(tmpdata.label, ft_channelselection('EOG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.eogscale;
  end
end
if ~isempty(cfg.ecgscale)
  scaled = 1;
  chansel = match_str(tmpdata.label, ft_channelselection('ECG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.ecgscale;
  end
end
if ~isempty(cfg.emgscale)
  scaled = 1;
  chansel = match_str(tmpdata.label, ft_channelselection('EMG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.emgscale;
  end
end
if ~isempty(cfg.megscale)
  scaled = 1;
  chansel = match_str(tmpdata.label, ft_channelselection('MEG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.megscale;
  end
end
if ~isempty(cfg.gradscale)
  scaled = 1;
  chansel = match_str(tmpdata.label, ft_channelselection('MEGGRAD', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.gradscale;
  end
end
if ~isempty(cfg.magscale)
  scaled = 1;
  chansel = match_str(tmpdata.label, ft_channelselection('MEGMAG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.magscale;
  end
end

if strcmp(cfg.method, 'channel')
  if scaled
    fprintf('showing the scaled data per channel, all trials at once\n');
  else
    fprintf('showing the data per channel, all trials at once\n');
  end
  [chansel, trlsel, cfg] = rejectvisual_channel(cfg, tmpdata);
elseif strcmp(cfg.method, 'trial')
  if scaled
    fprintf('showing the scaled per trial, all channels at once\n');
  else
    fprintf('showing the data per trial, all channels at once\n');
  end
  [chansel, trlsel, cfg] = rejectvisual_trial(cfg, tmpdata);
elseif strcmp(cfg.method, 'summary')
  if scaled
    fprintf('showing a summary of the scaled data for all channels and trials\n');
  else
    fprintf('showing a summary of the data for all channels and trials\n');
  end
  
  [chansel, trlsel, cfg] = rejectvisual_summary(cfg, tmpdata);
end

fprintf('%d trials marked as GOOD, %d trials marked as BAD\n', sum(trlsel), sum(~trlsel));
fprintf('%d channels marked as GOOD, %d channels marked as BAD\n', sum(chansel), sum(~chansel));

% construct an artifact matrix from the trl matrix
if isfield(data, 'sampleinfo')
  cfg.artifact = data.sampleinfo(~trlsel,:);
else
  % since sample numbers are unknown, it is not possible to remember them here
  cfg.artifact = [];
end

% remove artifacts from trl-matrix if present (but do *not* reconstruct the trl)
if isfield(cfg, 'trl')
  cfg.trl = cfg.trl(trlsel,:);
end

% show the user which trials are removed
removed = find(~trlsel);
if ~isempty(removed)
  fprintf('the following trials were removed: ');
  for i=1:(length(removed)-1)
    fprintf('%d, ', removed(i));
  end
  fprintf('%d\n', removed(end));
else
  fprintf('no trials were removed\n');
end

% remove the selected trials from the data
data.time  = data.time(trlsel);
data.trial = data.trial(trlsel);
if isfield(data, 'trialinfo'), data.trialinfo = data.trialinfo(trlsel,:); end;
if isfield(data, 'sampleinfo'),  data.sampleinfo  = data.sampleinfo(trlsel,:);  end;

% remove the offset vector if present (only applies to datasets that have been preprocessed a long time ago)
if isfield(data, 'offset')
  data = rmfield(data, 'offset');
end

if ~all(chansel)
  switch cfg.keepchannel
    case 'no'
      % show the user which channels are removed
      removed = find(~chansel);
      fprintf('the following channels were removed: ');
      for i=1:(length(removed)-1)
        fprintf('%s, ', data.label{removed(i)});
      end
      fprintf('%s\n', data.label{removed(end)});
      
      % remove channels that are not selected
      for i=1:length(data.trial)
        data.trial{i} = data.trial{i}(chansel,:);
      end
      data.label = data.label(chansel);
    case 'nan'
      % show the user which channels are removed
      removed = find(~chansel);
      fprintf('the following channels were filled with NANs: ');
      for i=1:(length(removed)-1)
        fprintf('%s, ', data.label{removed(i)});
      end
      fprintf('%s\n', data.label{removed(end)});
      
      % fill the data from the bad channels with nans
      for i=1:length(data.trial)
        data.trial{i}(~chansel,:) = nan;
      end
    case 'yes'
      % keep all channels, also when they are not selected
  end
end

% convert back to input type if necessary
switch dtype
  case 'timelock'
    data = ft_checkdata(data, 'datatype', 'timelock');
  otherwise
    % keep the output as it is
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous data
ft_postamble history data
ft_postamble savevar data
