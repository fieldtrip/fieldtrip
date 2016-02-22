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
%   cfg.method      = string, describes how the data should be shown, this can be
%                     'summary'  show a single number for each channel and trial (default)
%                     'channel'  show the data per channel, all trials at once
%                     'trial'    show the data per trial, all channels at once
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'),
%                     see FT_CHANNELSELECTION for details
%   cfg.keepchannel = string, determines how to deal with channels that are not selected, can be
%                     'no'          completely remove deselected channels from the data (default)
%                     'yes'         keep deselected channels in the output data
%                     'nan'         fill the channels that are deselected with NaNs
%                     'repair'      repair the deselected channels using FT_CHANNELREPAIR
%   cfg.neighbours  = neighbourhood structure, see also FT_PREPARE_NEIGHBOURS (required for repairing channels)
%   cfg.trials      = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.keeptrial   = string, determines how to deal with trials that are
%                     not selected, can be
%                     'no'     completely remove deselected trials from the data (default)
%                     'yes'    keep deselected trials in the output data
%                     'nan'    fill the trials that are deselected with NaNs
%   cfg.metric      = string, describes the metric that should be computed in summary mode
%                     for each channel in each trial, can be
%                     'var'       variance within each channel (default)
%                     'min'       minimum value in each channel
%                     'max'       maximum value each channel
%                     'maxabs'    maximum absolute value in each channel
%                     'range'     range from min to max in each channel
%                     'kurtosis'  kurtosis, i.e. measure of peakedness of the amplitude distribution
%                     'zvalue'    mean and std computed over all time and trials, per channel
%   cfg.latency     = [begin end] in seconds, or 'minperlength', 'maxperlength',
%                     'prestim', 'poststim' (default = 'maxperlength')
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
% To facilitate data-handling and distributed computing you can use
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
% Copyright (C) 2006-2016, Robert Oostenveld
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

% Undocumented options
% cfg.plotlayout = 'square' (default) or '1col', plotting every channel/trial under each other
% cfg.viewmode   = 'remove', 'toggle' or 'hide', only applies to summary mode (default = 'remove')

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% store original datatype
dtype = ft_datatype(data);

% check if the input data is valid for this function, this will convert it to raw if needed
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval',  {'metric',  'absmax',  'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method',  'absmax',  'maxabs'});

% resolve some common typing errors
cfg = ft_checkconfig(cfg, 'renamed',  {'keeptrials',  'keeptrial'});
cfg = ft_checkconfig(cfg, 'renamed',  {'keepchannels',  'keepchannel'});

% set the defaults
cfg.channel     = ft_getopt(cfg, 'channel'    , 'all');
cfg.trials      = ft_getopt(cfg, 'trials'     , 'all', 1);
cfg.latency     = ft_getopt(cfg, 'latency'    , 'maxperlength');
cfg.keepchannel = ft_getopt(cfg, 'keepchannel', 'no');
cfg.keeptrial   = ft_getopt(cfg, 'keeptrial'  , 'no');
cfg.feedback    = ft_getopt(cfg, 'feedback'   , 'textbar');
cfg.method      = ft_getopt(cfg, 'method'     , 'summary');
cfg.metric      = ft_getopt(cfg, 'metric'     , 'var');
cfg.alim        = ft_getopt(cfg, 'alim'       );
cfg.eegscale    = ft_getopt(cfg, 'eegscale'   );
cfg.eogscale    = ft_getopt(cfg, 'eogscale'   );
cfg.ecgscale    = ft_getopt(cfg, 'ecgscale'   );
cfg.emgscale    = ft_getopt(cfg, 'emgscale'   );
cfg.megscale    = ft_getopt(cfg, 'megscale'   );
cfg.gradscale   = ft_getopt(cfg, 'gradscale'  );
cfg.magscale    = ft_getopt(cfg, 'magscale'   );
cfg.plotlayout  = ft_getopt(cfg, 'plotlayout' , 'square');
cfg.viewmode    = ft_getopt(cfg, 'viewmode'   , 'remove');

% ensure that the preproc specific options are located in the cfg.preproc substructure
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});

% check required fields at the start, rather than further down in the code
if strcmp(cfg.keepchannel, 'repair')
  cfg = ft_checkconfig(cfg, 'required', 'neighbours');
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

% apply scaling to the selected channel types to equate the absolute numbers (i.e. fT and uV)
% make a seperate copy to prevent the original data from being scaled
tmpdata = data;
scaled  = false;
if ~isempty(cfg.eegscale)
  scaled = true;
  chansel = match_str(tmpdata.label, ft_channelselection('EEG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.eegscale;
  end
end
if ~isempty(cfg.eogscale)
  scaled = true;
  chansel = match_str(tmpdata.label, ft_channelselection('EOG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.eogscale;
  end
end
if ~isempty(cfg.ecgscale)
  scaled = true;
  chansel = match_str(tmpdata.label, ft_channelselection('ECG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.ecgscale;
  end
end
if ~isempty(cfg.emgscale)
  scaled = true;
  chansel = match_str(tmpdata.label, ft_channelselection('EMG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.emgscale;
  end
end
if ~isempty(cfg.megscale)
  scaled = true;
  chansel = match_str(tmpdata.label, ft_channelselection('MEG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.megscale;
  end
end
if ~isempty(cfg.gradscale)
  scaled = true;
  chansel = match_str(tmpdata.label, ft_channelselection('MEGGRAD', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.gradscale;
  end
end
if ~isempty(cfg.magscale)
  scaled = true;
  chansel = match_str(tmpdata.label, ft_channelselection('MEGMAG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.magscale;
  end
end

switch cfg.method
  case 'channel'
    if scaled
      fprintf('showing the scaled data per channel, all trials at once\n');
    else
      fprintf('showing the data per channel, all trials at once\n');
    end
    [chansel, trlsel, cfg] = rejectvisual_channel(cfg, tmpdata);

  case 'trial'
    if scaled
      fprintf('showing the scaled per trial, all channels at once\n');
    else
      fprintf('showing the data per trial, all channels at once\n');
    end
    [chansel, trlsel, cfg] = rejectvisual_trial(cfg, tmpdata);

  case 'summary'
    if scaled
      fprintf('showing a summary of the scaled data for all channels and trials\n');
    else
      fprintf('showing a summary of the data for all channels and trials\n');
    end
    [chansel, trlsel, cfg] = rejectvisual_summary(cfg, tmpdata);

  otherwise
    error('unsupported method %s', cfg.method);
end % switch method

fprintf('%d trials marked as GOOD, %d trials marked as BAD\n', sum(trlsel), sum(~trlsel));
fprintf('%d channels marked as GOOD, %d channels marked as BAD\n', sum(chansel), sum(~chansel));


if ~all(chansel)
  switch cfg.keepchannel
    case 'yes'
      % keep all channels, also when they are not selected
      fprintf('no channels were removed from the data\n');

    case 'no'
      % show the user which channels are removed
      removed = find(~chansel);
      fprintf('the following channels were removed: ');
      for i=1:(length(removed)-1)
        fprintf('%s, ', data.label{removed(i)});
      end
      fprintf('%s\n', data.label{removed(end)});

    case 'nan'
      % show the user which channels are removed
      removed = find(~chansel);
      fprintf('the following channels were filled with NANs: ');
      for i=1:(length(removed)-1)
        fprintf('%s, ', data.label{removed(i)});
      end
      fprintf('%s\n', data.label{removed(end)});
      % mark the selection as nan
      for i=1:length(data.trial)
        data.trial{i}(~chansel,:) = nan;
      end

    case 'repair'
      % show which channels are to be repaired
      removed = find(~chansel);
      fprintf('the following channels were repaired using FT_CHANNELREPAIR: ');
      for i=1:(length(removed)-1)
        fprintf('%s, ', data.label{removed(i)});
      end
      fprintf('%s\n', data.label{removed(end)});

      % create cfg struct for call to ft_channelrepair
      orgcfg = cfg;
      tmpcfg = [];
      tmpcfg.trials = 'all';
      tmpcfg.badchannel = data.label(~chansel);
      tmpcfg.neighbours = cfg.neighbours;
      if isfield(cfg, 'grad')
          tmpcfg.grad = cfg.grad;
      end
      if isfield(cfg, 'elec')
          tmpcfg.elec = cfg.elec;
      end
      % repair the channels that were selected as bad
      data = ft_channelrepair(tmpcfg, data);
      % restore the provenance information
      [cfg, data] = rollback_provenance(cfg, data);
      % restore the original trials parameter, it should not be 'all'
      cfg = copyfields(orgcfg, cfg, {'trials'});

    otherwise
      error('invalid specification of cfg.keepchannel')
  end % case
end % if ~all(chansel)

if ~all(trlsel)
  switch cfg.keeptrial
    case 'yes'
      % keep all trials, also when they are not selected
      fprintf('no trials were removed from the data\n');

    case 'no'
      % show the user which channels are removed
      removed = find(~trlsel);
      fprintf('the following trials were removed: ');
      for i=1:(length(removed)-1)
        fprintf('%d, ', removed(i));
      end
      fprintf('%d\n', removed(end));

    case 'nan'
      % show the user which trials are removed
      removed = find(~trlsel);
      fprintf('the following trials were filled with NANs: ');
      for i=1:(length(removed)-1)
        fprintf('%d, ', removed(i));
      end
      fprintf('%d\n', removed(end));
      % mark the selection as nan
      for i=removed
        data.trial{i}(:,:) = nan;
      end

    otherwise
      error('invalid specification of cfg.keeptrial')
  end % case
end % if ~all(trlsel)

if isfield(data, 'sampleinfo')
  % construct the matrix with sample numbers prior to making the selection
  cfg.artfctdef.(cfg.method).artifact = data.sampleinfo(~trlsel,:);
end

% perform the selection of channels and trials
orgcfg = cfg;
tmpcfg = [];
if strcmp(cfg.keepchannel, 'no')
  tmpcfg.channel = find(chansel);
end
if strcmp(cfg.keeptrial, 'no')
  tmpcfg.trials = find(trlsel); % note that it is keeptrial without S and trials with S
end
data = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);
% restore the original channels and trials parameters
cfg = copyfields(orgcfg, cfg, {'channel', 'trials'});

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
ft_postamble previous   data
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
