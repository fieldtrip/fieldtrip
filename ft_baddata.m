function [data] = ft_baddata(cfg, data)

% FT_BADDATA identifies bad data in a MEG or EEG dataset by looping over all trials
% and all channels. Each channel in each trial is considered separately, in the
% remainder of the help we will refer to this as "traces". Different methods are
% implemented, these are largely shared with those implemented in FT_REJECTVISUAL
% with the "summary" method. The methods are shortly described in detail below. Bad
% traces are replaced in the output data with nan.
%
% VAR, STD, MIN, MAX, MAXABS, RANGE, KURTOSIS, ZVALUE - compute the specified metric
% for each channel in each trial and check whether it exceeds the threshold.
%
% NEIGHBEXPVAR - identifies channels that cannot be explained very well by a linear
% combination of their neighbours. A general linear model is used to compute the
% explained variance. A value close to 1 means that a channel is similar to its
% neighbours, a value close to 0 indicates a "bad" channel.
%
% Use as
%   [data_clean] = ft_baddata(cfg, data)
% where the input data corresponds to the output from FT_PREPROCESSING.
%
% The configuration should contain
%   cfg.metric        = string, describes the metric that should be computed in summary mode for each channel in each trial, can be
%                       'var'          variance within each channel  (default)
%                       'std'          standard deviation within each channel 
%                       'db'           decibel value within each channel
%                       'mad'          median absolute deviation within each channel
%                       '1/var'        inverse variance within each channel
%                       'min'          minimum value in each channel
%                       'max'          maximum value in each channel
%                       'maxabs'       maximum absolute value in each channel
%                       'range'        range from min to max in each channel
%                       'kurtosis'     kurtosis, i.e. measure of peakedness of the amplitude distribution in trace
%                       'zvalue'       mean and std computed over all time and trials, per channel
%                       'neighbexpvar' relative variance explained by neighboring channels in each trial
%   cfg.threshold     = scalar, the appropriate value depends on the data characteristics and the metric
%   cfg.feedback      = 'yes' or 'no', whether to show an image of the neighbour values (default = 'no')
%
% The following options allow you to make a pre-selection
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.trials      = 'all' or a selection given as a 1xN vector (default = 'all')
%
% See also FT_BADCHANNEL, FT_BADSEGMENT, FT_REJECTVISUAL, FT_CHANNELREPAIR

% Undocumented options
%   cfg.thresholdside = above or below

% Copyright (C) 2024, Robert Oostenveld
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
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels', 'trial'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'required',   'metric');

% ensure that the preproc specific options are located in the cfg.preproc substructure
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});

% set the defaults
cfg.channel       = ft_getopt(cfg, 'channel', 'all');
cfg.trials        = ft_getopt(cfg, 'trials', 'all', true);
cfg.neighbours    = ft_getopt(cfg, 'neighbours');
cfg.nbdetect      = ft_getopt(cfg, 'nbdetect', 'median');
cfg.feedback      = ft_getopt(cfg, 'feedback', 'no');
cfg.thresholdside = ft_getopt(cfg, 'thresholdside', []); % the default depends on cfg.metric, see below

if isempty(cfg.thresholdside)
  if ismember(cfg.metric, {'var', 'std', 'db', 'mad', '1/var', 'max', 'maxabs', 'range', 'kurtosis', 'zvalue', 'maxzvalue', 'neighbstdratio'})
    % large positive values indicate an artifact, so check for values ABOVE the threshold
    cfg.thresholdside = 'above';
  elseif ismember(cfg.metric, {'min', 'neighbexpvar', 'neighbcorr'})
    % very negative values or small positive values indicate an artifact, so check for values BELOW the threshold
    cfg.thresholdside = 'below';
  else
    % there are also a few where one could look at either side, these require the user to make a choice
    ft_error('you must specify cfg.thresholdside');
  end
end

% select trials and channels of interest
tmpcfg = keepfields(cfg, {'trials', 'channel', 'tolerance', 'latency', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
data   = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

ntrl  = length(data.trial);
nchan = length(data.label);

if contains(cfg.metric, 'zvalue')
  % cellmean and cellstd (see FT_DENOISE_PCA) would work instead of for-loops, but they are too memory-intensive
  runsum = zeros(nchan, 1);
  runss  = zeros(nchan, 1);
  runnum = 0;
  for trl=1:ntrl
    dat = preproc(data.trial{trl}, data.label, data.time{trl}, cfg.preproc);
    runsum = runsum + nansum(dat, 2);
    runss  = runss  + nansum(dat.^2, 2);
    runnum = runnum + sum(isfinite(dat), 2);
  end
  mval = runsum./runnum;
  sd   = sqrt(runss./runnum - (runsum./runnum).^2);
else
  mval = [];
  sd   = [];
end

if contains(cfg.metric, 'neighb')
  cfg = ft_checkconfig(cfg, 'required', 'neighbours');
  % creates a NxN Boolean matrix that describes whether channels are connected as neighbours
  connectivity = channelconnectivity(cfg, data);
else
  connectivity = [];
end

% compute the artifact value for each trial and each channel
level = nan(nchan, ntrl);
for trl=1:ntrl
  dat = preproc(data.trial{trl}, data.label, data.time{trl}, cfg.preproc);
  level(:,trl) = artifact_level(dat, cfg.metric, mval, sd, connectivity);
end

% find channels and trials with a value that exceeds the threshold
switch cfg.thresholdside
  case 'below'
    bad = level<cfg.threshold;
  case 'above'
    bad = level>cfg.threshold;
end

ft_info('identified %d out of %d traces as bad (%.0f %%)\n', sum(bad(:)), length(bad(:)), 100*mean(bad(:)));

for trl=1:ntrl
  for chan=1:nchan
    if bad(chan,trl)
      data.trial{trl}(chan,:) = nan;
    end
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   data
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
