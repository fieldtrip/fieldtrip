function [cfg] = ft_badchannel(cfg, data)

% FT_BADCHANNEL tries to identify bad channels in a MEG or EEG dataset. Different
% methods are implemented to identify bad channels, these are largely shared with
% those implemented in FT_REJECTVISUAL with the summary method. The methods are
% shortly described in detail below.
%
% VAR, STD, MIN, MAX, MAXABS, RANGE, KURTOSIS, ZVALUE - compute the specified metric
% for each channel in each trial and check whether it exceeds the threshold.
%
% NEIGHBEXPVAR - identifies channels that cannot be explained very well by a linear
% combination of their neighbours. A general linear model is used to compute the
% explained variance. A value close to 1 means that a channel is similar to its
% neighbours, a value close to 0 indicates a "bad" channel.
%
% NEIGHBCORR - identifies channels that have low correlation with each of their
% neighbours. The rationale is that "bad" channel have inherent noise that is
% uncorrelated with other sensors.
%
% NEIGHBSTDRATIO - identifies channels that have a standard deviation which is very
% different from that of each of their neighbours. This computes the difference in
% the standard deviation of each channel to each of its neighbours, relative to that
% of the neighbours.
%
% Use as
%   [cfg] = ft_badchannel(cfg, data)
% where the input data corresponds to the output from FT_PREPROCESSING.
%
% The configuration should contain
%   cfg.metric        = string, describes the metric that should be computed in summary mode for each channel in each trial, can be
%                       'var'       variance within each channel (default)
%                       'min'       minimum value in each channel
%                       'max'       maximum value in each channel
%                       'maxabs'    maximum absolute value in each channel
%                       'range'     range from min to max in each channel
%                       'kurtosis'  kurtosis, i.e. measure of peakedness of the amplitude distribution
%                       'zvalue'    mean and std computed over all time and trials, per channel
%   cfg.threshold     = scalar, the optimal value depends on the methods and on the data characteristics
%   cfg.neighbours    = neighbourhood structure, see FT_PREPARE_NEIGHBOURS for details
%   cfg.nbdetect      = 'any', 'most', 'all', 'median', see below (default = 'median')
%   cfg.feedback      = 'yes' or 'no', whether to show an image of the neighbour values (default = 'no')
%
% The following options allow you to make a pre-selection
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.trials      = 'all' or a selection given as a 1xN vector (default = 'all')
%
% The 'neighcorrel' and 'neighstdratio' methods implement the bad channel detection
% (more or less) according to the paper "Adding dynamics to the Human Connectome
% Project with MEG", Larson-Prior et al. https://doi.org/10.1016/j.neuroimage.2013.05.056.
%
% Most methods compute a scalar value for each channel that can simply be
% thresholded. The NEIGHBCORR and NEIGHBSTDRATIO compute a vector with a value for
% each of the neighbour of a channel. The cfg.nbdetect option allows you to specify
% whether you want to flag the channel as bad in case 'all' of its neighbours exceed
% the threshold, if 'most' exceed the threshold, or if 'any' of them exceeds the
% threshold. Note that when you specify 'any', then all channels neighbouring a bad
% channel will also be marked as bad, since they all have at least one bad neighbour.
% You can also specify 'median', in which case the threshold is applied to the median
% value over neighbours.
%
% See also FT_BADSEGMENT, FT_REJECTVISUAL, FT_CHANNELREPAIR

% Undocumented options
%   cfg.thresholdside = above or below

% Copyright (C) 2021, Robert Oostenveld
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
  if ismember(cfg.metric, {'var', 'std', 'max', 'maxabs', 'range', 'kurtosis', '1/var', 'zvalue', 'maxzvalue', 'neighbstdratio'})
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
tmpcfg = keepfields(cfg, {'trials', 'channel', 'tolerance', 'latency', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
data   = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

ntrl  = length(data.trial);
nchan = length(data.label);
badchannel = false(nchan,1);

if contains(cfg.metric, 'zvalue')
  % cellmean and cellstd (see FT_DENOISE_PCA) would work instead of for-loops, but they are too memory-intensive
  runsum = zeros(nchan, 1);
  runss  = zeros(nchan, 1);
  runnum = 0;
  for chan=1:ntrl
    dat = preproc(data.trial{chan}, data.label, data.time{chan}, cfg.preproc);
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

for trl=1:ntrl
  % compute the artifact value for each channel in this trial
  level = artifact_level(data.trial{trl}, cfg.metric, mval, sd, connectivity);
  
  if isvector(level)
    % find channels with a value that exceeds the threshold
    switch cfg.thresholdside
      case 'below'
        badchannel = badchannel | level<cfg.threshold;
      case 'above'
        badchannel = badchannel | level>cfg.threshold;
    end
    
  else
    % identify channels with one of their neighbours values that exceeds the threshold
    for chan=1:nchan
      nblevel = level(chan,:);            % select this channel from the matrix
      nblevel = nblevel(~isnan(nblevel)); % only select its actual neighbours
      switch cfg.nbdetect
        case 'all'
          switch cfg.thresholdside
            case 'below'
              badchannel(chan) = badchannel(chan) | all(nblevel<cfg.threshold);
            case 'above'
              badchannel(chan) = badchannel(chan) | all(nblevel>cfg.threshold);
          end % switch
        case 'most'
          switch cfg.thresholdside
            case 'below'
              badchannel(chan) = badchannel(chan) | most(nblevel<cfg.threshold);
            case 'above'
              badchannel(chan) = badchannel(chan) | most(nblevel>cfg.threshold);
          end % switch
        case 'any'
          switch cfg.thresholdside
            case 'below'
              badchannel(chan) = badchannel(chan) | any(nblevel<cfg.threshold);
            case 'above'
              badchannel(chan) = badchannel(chan) | any(nblevel>cfg.threshold);
          end % switch
        case 'median'
          switch cfg.thresholdside
            case 'below'
              badchannel(chan) = badchannel(chan) | nanmedian(nblevel,2)<cfg.threshold;
            case 'above'
              badchannel(chan) = badchannel(chan) | nanmedian(nblevel,2)>cfg.threshold;
          end % switch
        otherwise
          ft_error('incorrect specification of cfg.nbdetect');
      end
    end % for each channel
    
  end % if isvector
end % for each trial

ft_info('identified %d out of %d channels as bad\n', sum(badchannel), length(badchannel));

% keep track of bad channels
cfg.badchannel = data.label(badchannel);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous data
ft_postamble provenance
ft_postamble hastoolbox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = most(x)
tf = sum(x(:)==true)>(numel(x)/2);
