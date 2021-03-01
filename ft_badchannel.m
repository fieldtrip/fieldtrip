function cfg = ft_badchannel(cfg, data)

% FT_BADCHANNEL tries to identify bad channels in a MEG or EEG dataset. Different
% methods are implemented to identify bad channels, these are largely shared with the
% methods implemented in FT_REJECTVISUAL with the summary method. The methods are
% shortly described in detail below.
%
% VAR, STD, MIN, MAX, RANGE, KURTOSIS - concatenate the data over trials, compute the
% specified metric for each channel and check whether it exceeds the threshold.
%
% NEIGHBCORR - identifies channels that have low correlation with their neighbours.
% The rationale is that "bad" channel have inherent noise that is uncorrelated with
% other sensors. The average correlation of MEG magnetometer channels in a BTI248 MEG
% system with their neighbours within a distance of 4 cm is between 0.6 and 0.8.
%
% NEIGHBSTDRATIO - identifies channels that have a standard deviation which is very
% different from their neighbours. This computes the difference in the standard
% deviation of each channel to each of its neighbours, relative to that of the
% neighbours. A value of 0.5 has empirically been determined as adequate for
% identifying bad channels in a BTI248 MEG system.
%
% Use as
%   cfg = ft_badchannel(cfg, data)
% where the input data corresponds to the output from FT_PREPROCESSING.
%
% The configuration should contain
%   cfg.method        = 'var', 'std', 'min', 'max', 'range', 'kurtosis', 'neighbcorr', 'neighbstdratio'
%   cfg.threshold     = scalar
%   cfg.thresholdside = 'above' or 'below' (default is automatic)
%   cfg.neighbours    = neighbourhood structure, see FT_PREPARE_NEIGHBOURS for details
%   cfg.nbdetect      = 'any', 'most', 'all', 'median', see below (default = 'all')
%   cfg.feedback      = 'yes' or 'no', whether to show an image of the neighbour values (default = 'no')
%
% The 'neighcorrel' and 'neighstdratio' methods implement the bad channel detection
% according to the paper "Adding dynamics to the Human Connectome Project with MEG",
% Larson-Prior et al. https://doi.org/10.1016/j.neuroimage.2013.05.056
%
% The cfg.nbdetect option allows you to specify whether you want to flag a channel
% as bad in case 'all' of its neighbours exceed the threshold, if 'most' exceed the
% threshold, or if 'any' of them exceeds the threshold. When you specify 'any', then
% all channels neighbouring a bad channel will also be marked as bad, since they all
% have at least one bad neighbour. Alternatively, you can specify 'median', in which
% case the threshold is applied to the median value over neighbours.
%
% See also FT_CHANNELREPAIR, FT_REJECTVISUAL

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
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

% ensure that the configuration is consistent
cfg = ft_checkconfig(cfg, 'required', 'method');

% set the defaults
cfg.nbdetect      = ft_getopt(cfg, 'nbdetect', 'all');
cfg.feedback      = ft_getopt(cfg, 'feedback', 'no');
cfg.thresholdside = ft_getopt(cfg, 'thresholdside', []); % the default depends on cfg.method, see below

if isempty(cfg.thresholdside)
  if ismember(cfg.method, {'var', 'std', 'max', 'maxabs', 'range', 'kurtosis', '1/var'})
    % large positive values indicate an artifact, so check for values ABOVE the threshold
    cfg.thresholdside = 'above';
  elseif ismember(cfg.method, {'neighbexpvar', 'neighbcorr'})
    % small positive values indicate an artifact, so check for values BELOW the threshold
    cfg.thresholdside = 'below';
  elseif ismember(cfg.method, {'min'})
    % large negative values indicate an artifact, so check for values BELOW the threshold
    cfg.thresholdside = 'below';
  else
    % there are also a few where ohne could look at either side, these require the user to make a choice
    ft_error('cfg.thresholdside is required');
  end
end

if startsWith(cfg.method, 'neighb')
  cfg = ft_checkconfig(cfg, 'required', 'neighbours');
  % creates a NxN Boolean matrix that describes whether channels are connected as neighbours
  connectivity = channelconnectivity(cfg, data);
else
  connectivity = [];
end

% select channels and trials of interest
tmpcfg = keepfields(cfg, {'channel', 'trials', 'tolerance', 'showcallinfo'});
data = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

ft_info('concatenating data');
tim = cat(2, data.time{:});
dat = cat(2, data.trial{:});
lab = data.label;
nchan = length(lab);

% compute the artifact values
level = artifact_level(dat, cfg.method, [], [], connectivity);

if isvector(level)
  % give some feedback
  ft_info('min(range)    = %g\n', min(level));
  ft_info('max(range)    = %g\n', max(level));
  ft_info('mean(range)   = %g\n', mean(level));
  ft_info('median(range) = %g\n', median(level));
  ft_info('mad(range)    = %g\n', mad(level));
  
  % find channels with a level below or above the threshold
  switch cfg.thresholdside
    case 'below'
      bad = level<cfg.threshold;
    case 'above'
      bad = level>cfg.threshold;
  end
  
elseif ismatrix(level)
  % give some feedback
  tmp = level(~isnan(level(:)));
  ft_info('min(range)    = %g\n', min(tmp));
  ft_info('max(range)    = %g\n', max(tmp));
  ft_info('mean(range)   = %g\n', mean(tmp));
  ft_info('median(range) = %g\n', median(tmp));
  ft_info('mad(range)    = %g\n', mad(tmp));
  if istrue(cfg.feedback)
    figure; histogram(tmp, round(length(tmp)/10));
    xy = axis;
    line([cfg.threshold cfg.threshold], xy(3:4));
    figure; imagesc(level); colorbar
  end
  
  % identify channels with a neighbour value ABOVE the threshold
  switch cfg.nbdetect
    case 'all'
      bad = false(1,nchan);
      for i=1:nchan
        tmp = level(i,:);
        tmp = tmp(~isnan(tmp)); % only select the actual neighbours
        switch cfg.thresholdside
          case 'below'
            bad(i) = all(tmp<cfg.threshold);
          case 'above'
            bad(i) = all(tmp>cfg.threshold);
        end % switch
      end
    case 'most'
      bad = false(1,nchan);
      for i=1:nchan
        tmp = level(i,:);
        tmp = tmp(~isnan(tmp)); % only select the actual neighbours
        switch cfg.thresholdside
          case 'below'
            bad(i) = most(tmp<cfg.threshold);
          case 'above'
            bad(i) = most(tmp>cfg.threshold);
        end % switch
      end
    case 'any'
      bad = false(1,nchan);
      for i=1:nchan
        tmp = level(i,:);
        tmp = tmp(~isnan(tmp)); % only select the actual neighbours
        switch cfg.thresholdside
          case 'below'
            bad(i) = any(tmp<cfg.threshold);
          case 'above'
            bad(i) = any(tmp>cfg.threshold);
        end % switch
      end
    case 'median'
      switch cfg.thresholdside
        case 'below'
          bad = nanmedian(level,2)<cfg.threshold;
        case 'above'
          bad = nanmedian(level,2)>cfg.threshold;
      end % switch
    otherwise
      ft_error('incorrect specification of cfg.nbdetect');
  end
  
end % if isvector or ismatrix

ft_info('identified %d bad channels from %d channels\n', sum(bad), length(bad));
cfg.badchannel = lab(bad);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data
ft_postamble provenance
ft_postamble hastoolbox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = most(x)
tf = sum(x(:)==true)>(numel(x)/2);
