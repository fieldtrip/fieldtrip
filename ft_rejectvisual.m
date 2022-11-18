function [data] = ft_rejectvisual(cfg, data)

% FT_REJECTVISUAL shows the preprocessed data in all channels and/or trials to allow
% the user to make a visual selection of the data that should be rejected. The data
% can be displayed in a "summary" mode, in which case the variance (or another
% metric) in each channel and each trial is computed. Alternatively, all channels can
% be shown at once allowing paging through the trials, or all trials can be shown,
% allowing paging through the channels.
%
% Use as
%   [data] = ft_rejectvisual(cfg, data)
%
% The configuration can contain
%   cfg.method      = string, describes how the data should be shown, this can be
%                     'summary'  show a single number for each channel and trial (default)
%                     'channel'  show the data per channel, all trials at once
%                     'trial'    show the data per trial, all channels at once
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.keepchannel = string, determines how to deal with channels that are not selected, can be
%                     'no'          completely remove deselected channels from the data (default)
%                     'yes'         keep deselected channels in the output data
%                     'nan'         fill the channels that are deselected with NaNs
%                     'zero'        fill the channels that are deselected with zeros
%                     'repair'      repair the deselected channels using FT_CHANNELREPAIR
%   cfg.trials      = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.keeptrial   = string, determines how to deal with trials that are
%                     not selected, can be
%                     'no'     completely remove deselected trials from the data (default)
%                     'yes'    keep deselected trials in the output data
%                     'nan'    fill the trials that are deselected with NaNs
%                     'zero'   fill the trials that are deselected with zeros
%   cfg.metric      = string, describes the metric that should be computed in summary mode
%                     for each channel in each trial, can be
%                     'var'       variance within each channel (default)
%                     'min'       minimum value in each channel
%                     'max'       maximum value each channel
%                     'maxabs'    maximum absolute value in each channel
%                     'range'     range from min to max in each channel
%                     'kurtosis'  kurtosis, i.e. measure of peakedness of the amplitude distribution
%                     'zvalue'    mean and std computed over all time and trials, per channel
%   cfg.neighbours  = neighbourhood structure, see FT_PREPARE_NEIGHBOURS for details
%   cfg.latency     = [begin end] in seconds, or 'all', 'minperiod', 'maxperiod', 'prestim', 'poststim' (default = 'all')
%   cfg.viewmode    = 'remove', 'toggle' or 'hide', only applies to summary mode (default = 'remove')
%   cfg.box         = string, 'yes' or 'no' whether to draw a box around each graph (default = 'no')
%   cfg.ylim        = 'maxmin', 'maxabs', 'zeromax', 'minzero', or [ymin ymax] (default = 'maxmin')
%
% The following options for the scaling of the EEG, EOG, ECG, EMG, MEG and NIRS channels
% is optional and can be used to bring the absolute numbers of the different
% channel types in the same range (e.g. fT and uV). The channel types are determined
% from the input data using FT_CHANNELSELECTION.
%   cfg.eegscale    = number, scaling to apply to the EEG channels prior to display
%   cfg.eogscale    = number, scaling to apply to the EOG channels prior to display
%   cfg.ecgscale    = number, scaling to apply to the ECG channels prior to display
%   cfg.emgscale    = number, scaling to apply to the EMG channels prior to display
%   cfg.megscale    = number, scaling to apply to the MEG channels prior to display
%   cfg.gradscale   = number, scaling to apply to the MEG gradiometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.magscale    = number, scaling to apply to the MEG magnetometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.nirsscale   = number, scaling to apply to the NIRS channels prior to display
%   cfg.mychanscale = number, scaling to apply to the channels specified in cfg.mychan
%   cfg.mychan      = Nx1 cell-array with selection of channels
%   cfg.chanscale   = Nx1 vector with scaling factors, one per channel specified in cfg.channel
%
% Optionally, the raw data is preprocessed (filtering etc.) prior to displaying it or
% prior to computing the summary metric. The preprocessing and the selection of the
% latency window is NOT applied to the output data.
%
% The following settings are useful for identifying EOG artifacts:
%   cfg.preproc.bpfilter    = 'yes'
%   cfg.preproc.bpfilttype  = 'but'
%   cfg.preproc.bpfreq      = [1 15]
%   cfg.preproc.bpfiltord   = 4
%   cfg.preproc.rectify     = 'yes'
%
% The following settings are useful for identifying muscle artifacts:
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

% Copyright (C) 2005-2006, Markus Bauer, Robert Oostenveld
% Copyright (C) 2006-2021, Robert Oostenveld
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

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% store original datatype
dtype = ft_datatype(data);

% check if the input data is valid for this function, this will convert it to raw if needed
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels', 'trial'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'renamedval', {'metric',  'absmax',  'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval', {'method',  'absmax',  'maxabs'});

% resolve some common typing errors
cfg = ft_checkconfig(cfg, 'renamed',  {'keeptrials',  'keeptrial'});
cfg = ft_checkconfig(cfg, 'renamed',  {'keepchannels',  'keepchannel'});
cfg = ft_checkconfig(cfg, 'renamed',  {'alim',  'ylim'});

% set the defaults
cfg.channel     = ft_getopt(cfg, 'channel'    , 'all');
cfg.trials      = ft_getopt(cfg, 'trials'     , 'all', true);
cfg.latency     = ft_getopt(cfg, 'latency'    , 'maxperiod');
cfg.keepchannel = ft_getopt(cfg, 'keepchannel', 'no');
cfg.keeptrial   = ft_getopt(cfg, 'keeptrial'  , 'no');
cfg.feedback    = ft_getopt(cfg, 'feedback'   , 'textbar');
cfg.method      = ft_getopt(cfg, 'method'     , 'summary');
cfg.metric      = ft_getopt(cfg, 'metric'     , 'var');
cfg.ylim        = ft_getopt(cfg, 'ylim'       );
cfg.eegscale    = ft_getopt(cfg, 'eegscale'   );
cfg.eogscale    = ft_getopt(cfg, 'eogscale'   );
cfg.ecgscale    = ft_getopt(cfg, 'ecgscale'   );
cfg.emgscale    = ft_getopt(cfg, 'emgscale'   );
cfg.megscale    = ft_getopt(cfg, 'megscale'   );
cfg.gradscale   = ft_getopt(cfg, 'gradscale'  );
cfg.magscale    = ft_getopt(cfg, 'magscale'   );
cfg.ylim        = ft_getopt(cfg, 'ylim'       , 'maxmin');
cfg.layout      = ft_getopt(cfg, 'layout'     , 'ordered');
cfg.viewmode    = ft_getopt(cfg, 'viewmode'   , 'remove');
cfg.box         = ft_getopt(cfg, 'box'        , 'no');

% this is needed for the figure title
if isfield(cfg, 'dataname') && ~isempty(cfg.dataname)
  cfg.dataname = cfg.dataname;
elseif isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
  cfg.dataname = cfg.inputfile;
elseif nargin>1
  cfg.dataname = inputname(2);
else
  cfg.dataname = {};
end

% ensure that the preproc specific options are located in the cfg.preproc substructure
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});

% check required fields at the start, rather than further down in the code
if strcmp(cfg.keepchannel, 'repair')
  cfg = ft_checkconfig(cfg, 'required', 'neighbours');
end

% apply scaling to the selected channel types to equate the absolute numbers (i.e. fT and uV)
fn = fieldnames(cfg);
tmpcfg = keepfields(cfg, fn(endsWith(fn, 'scale') | startsWith(fn, 'mychan') | strcmp(fn, 'channel')));
tmpcfg.parameter = 'trial';
tmpdata = chanscale_common(tmpcfg, data);
scaled = ~isequal(data.trial, tmpdata.trial);

% at this moment it is important that NO data selection is made, all data is passed through to the subfunctions
% which subsequently refine the initial cfg-based inclusion/exclusion of channels and trials
% (important because this way the original channel/trial indices are
% available in the GUI)

% to highlight to the user that cfg.trials/cfg.channel operate on the same
% selection of trials/channels as the user interface, mention here the
% consequences of the selection *before* any user interaction
ntrl_all = length(data.trial);
if isequal(cfg.trials, 'all') || isempty(cfg.trials)
  ntrl_keep = ntrl_all;
elseif isnumeric(cfg.trials)
  ntrl_keep = numel(cfg.trials);
elseif islogical(cfg.trials)
  ntrl_keep = sum(cfg.trials);
end
nchan_all  = numel(data.label);
nchan_keep = numel(ft_channelselection(cfg.channel, data.label));
fprintf('before GUI interaction: %d trials marked to INCLUDE, %d trials marked to EXCLUDE\n', ntrl_keep, ntrl_all-ntrl_keep);
fprintf('before GUI interaction: %d channels marked to INCLUDE, %d channels marked to EXCLUDE\n', nchan_keep, nchan_all-nchan_keep);

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
    ft_error('unsupported method %s', cfg.method);
end % switch method

fprintf('after GUI interaction: %d trials marked to INCLUDE, %d trials marked to EXCLUDE\n', sum(trlsel), sum(~trlsel));
fprintf('after GUI interaction: %d channels marked to INCLUDE, %d channels marked to EXCLUDE\n', sum(chansel), sum(~chansel));

% these are to be removed, filled with nan/zero, or kept in the output
badchannel = tmpdata.label(~chansel);
badsegment = find(~trlsel);

if ~all(chansel)
  switch cfg.keepchannel
    case 'yes'
      % keep all channels, also when they are not selected
      fprintf('no channels were removed from the data\n');

    case 'no'
      % show the user which channels are removed
      fprintf('the following channels were removed: ');

    case 'nan'
      % show the user which channels are nan-filled
      fprintf('the following channels were filled with NaNs: ');

      % mark the selection as nan
      for i=1:length(data.trial)
        data.trial{i}(~chansel,:) = nan;
      end

    case 'zero'
      % show the user which channels are zero-filled
      fprintf('the following channels were filled with zeros: ');

      % mark the selection as zero
      for i=1:length(data.trial)
        data.trial{i}(~chansel,:) = 0;
      end

    case 'repair'
      % create cfg struct for call to FT_CHANNELREPAIR
      orgcfg = cfg;
      tmpcfg = [];
      if isfield(data, 'grad') || isfield(data, 'elec') || isfield(data, 'opto')
        tmpcfg.method = 'weighted';
      else
        tmpcfg.method = 'average';
      end
      tmpcfg.trials = 'all'; % here we are only dealing with bad channels, bad trials will be removed further down (if applicable)
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
      [orgcfg, data] = rollback_provenance(orgcfg, data);
      % restore the original trials parameter, it should not be 'all'
      cfg = copyfields(orgcfg, cfg, {'trials'});

      % show which channels were repaired
      fprintf('the following channels were repaired using FT_CHANNELREPAIR: ');

    otherwise
      ft_error('invalid specification of cfg.keepchannel')
  end % case

  % provide the channel feedback
  if any(strcmp({'no', 'nan', 'repair'}, cfg.keepchannel))
    for i=1:(length(badchannel)-1)
      fprintf('%s, ', badchannel{i});
    end
    fprintf('%s\n', badchannel{end});
  end

end % if ~all(chansel)

if ~all(trlsel)
  switch cfg.keeptrial
    case 'yes'
      % keep all trials, also when they are not selected
      fprintf('no trials were removed from the data\n');

    case 'no'
      % show the user which channels are removed
      fprintf('the following trials were removed: ');

    case 'nan'
      % show the user which trials are nan-filled
      fprintf('the following trials were filled with NaNs: ');

      % mark the selection as nan
      for i = badsegment
        data.trial{i}(:,:) = nan;
      end

    case 'zero'
      % show the user which trials are zero-filled
      fprintf('the following trials were filled with zeros: ');

      % mark the selection as zero
      for i = badsegment
        data.trial{i}(:,:) = 0;
      end

    otherwise
      ft_error('invalid specification of cfg.keeptrial')
  end % case

  % provide the trial feedback
  if any(strcmp({'no', 'nan', 'repair'}, cfg.keeptrial))
    for i=1:(length(badsegment)-1)
      fprintf('%d, ', badsegment(i));
    end
    fprintf('%d\n', badsegment(end));
  end
end % if ~all(trlsel)

% keep track of bad segments
if isfield(data, 'sampleinfo')
  % this format is consistent with that of other artifact detection functions
  % construct the artifact matrix prior to making the selection
  cfg.artfctdef.(cfg.method).artifact = data.sampleinfo(badsegment,:);
end

% keep track of bad channels
cfg.badchannel = badchannel;

% these represent the channels and trials that are retained in the output
cfg.channel = data.label(chansel);
cfg.trials = find(trlsel);

% perform the actual selection of channels and trials
tmpcfg = [];
if strcmp(cfg.keepchannel, 'no')
  tmpcfg.channel = cfg.channel;
end
if strcmp(cfg.keeptrial, 'no')
  tmpcfg.trials = cfg.trials;
end
data = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

% convert back to input type if necessary
switch dtype
  case 'timelock'
    data = ft_checkdata(data, 'datatype', 'timelock');
  otherwise
    % keep the output as it is
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   data
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
