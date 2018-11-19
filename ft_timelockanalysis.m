function [timelock] = ft_timelockanalysis(cfg, data)

% FT_TIMELOCKANALYSIS computes the timelocked average ERP/ERF and
% optionally computes the covariance matrix. 
%
% Use as
%   [timelock] = ft_timelockanalysis(cfg, data)
%
% The data should be organised in a structure as obtained from the
% FT_PREPROCESSING function. The configuration should be according to
%
%   cfg.channel            = Nx1 cell-array with selection of channels (default = 'all'),
%                            see FT_CHANNELSELECTION for details
%   cfg.trials             = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.covariance         = 'no' or 'yes' (default = 'no')
%   cfg.covariancewindow   = 'prestim', 'poststim', 'all' or [begin end] (default = 'all')
%   cfg.keeptrials         = 'yes' or 'no', return individual trials or average (default = 'no')
%   cfg.removemean         = 'no' or 'yes' for covariance computation (default = 'yes')
%   cfg.vartrllength       = 0, 1 or 2 (see below)
%
% Depending on cfg.vartrllength, variable length trials and trials with
% differences in their time axes (so even if they are of the same length, e.g. 1
% second snippets of data cut from a single long recording) are treated differently:
%   0 - do not accept variable length trials [default]
%   1 - accept variable length trials, but only take those trials in which
%       data is present in both the average and the covariance window
%   2 - accept variable length trials, use all available trials
%       the available samples in every trial will be used for the
%       average and covariance computation. Missing values are replaced
%       by NaN and are not included in the computation.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_TIMELOCKGRANDAVERAGE, FT_TIMELOCKSTATISTICS

% Guidelines for use in an analysis pipeline:
% after FT_TIMELOCKANALYSIS you will have timelocked data - i.e., event-related
% fields (ERFs) or potentials (ERPs) - represented as the average and/or
% covariance over trials.
% This usually serves as input for one of the following functions:
%    * FT_TIMELOCKBASELINE      to perform baseline normalization
%    * FT_TIMELOCKGRANDAVERAGE  to compute the ERP/ERF average and variance over multiple subjects
%    * FT_TIMELOCKSTATISTICS    to perform parametric or non-parametric statistical tests
% Furthermore, the data can be visualised using the various plotting
% functions, including:
%    * FT_SINGLEPLOTER          to plot the ERP/ERF of a single channel or the average over multiple channels
%    * FT_TOPOPLOTER            to plot the topographic distribution over the head
%    * FT_MULTIPLOTER           to plot ERPs/ERFs in a topographical layout

% FIXME if input is one raw trial, the covariance is not computed correctly
%
% Undocumented local options:
% cfg.feedback
% cfg.preproc
%
% Deprecated options:
% cfg.latency
% cfg.blcovariance
% cfg.blcovariancewindow
% cfg.normalizevar
% cfg.normalizecov

% Copyright (C) 2018, Jan-Mathijs Schoffelen
% Copyright (C) 2003-2006, Markus Bauer
% Copyright (C) 2003-2006, Robert Oostenveld
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
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'normalizecov', 'normalizevar'});
cfg = ft_checkconfig(cfg, 'forbidden',  {'latency', 'blcovariance', 'blcovariancewindow'});
cfg = ft_checkconfig(cfg, 'renamed',     {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamed',     {'blcwindow', 'baselinewindow'});

% set the defaults
cfg.trials       = ft_getopt(cfg, 'trials',      'all', 1);
cfg.channel      = ft_getopt(cfg, 'channel',     'all');
cfg.keeptrials   = ft_getopt(cfg, 'keeptrials',  'no');
cfg.vartrllength = ft_getopt(cfg, 'vartrllength', 0);
cfg.feedback     = ft_getopt(cfg, 'feedback',     'text');
cfg.preproc      = ft_getopt(cfg, 'preproc',      []);
cfg.covariance       = ft_getopt(cfg, 'covariance',      'no');
cfg.covariancewindow = ft_getopt(cfg, 'covariancwindow', 'all');
cfg.removemean       = ft_getopt(cfg, 'removemean',      'yes');

% create logical flags for convenience
keeptrials = istrue(cfg.keeptrials);
computecov = istrue(cfg.covariance);

% ensure that the preproc specific options are located in the cfg.preproc substructure
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});

% select trials and channels of interest
tmpcfg = keepfields(cfg, {'trials', 'channel', 'latency'});
data   = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

if ~isempty(cfg.preproc)
  % preprocess the data, i.e. apply filtering, baselinecorrection, etc.
  fprintf('applying preprocessing options\n');
  if ~isfield(cfg.preproc, 'feedback')
    cfg.preproc.feedback = cfg.feedback;
  end
  data = ft_preprocessing(cfg.preproc, data);
  [cfg.preproc, data] = rollback_provenance(cfg.preproc, data);
end

% convert to a timelock structure with trials kept and NaNs for missing
% data points
data = ft_checkdata(data, 'datatype', 'timelock');

% check whether trial type (vartrllength) matches the user-specified type
if ~ismember(cfg.vartrllength, [0 1 2])
  ft_error('unknown value for vartrllength');
elseif cfg.vartrllength==0 && any(~isfinite(data.trial(:)))
  ft_error('data has variable trial lengths, you specified not to accept that');
elseif cfg.vartrllength==1 && all(isfinite(data.trial(:)))
  disp('data is of type fixed length !');
elseif cfg.vartrllength==2 && keeptrials
  disp('processing and keeping variable length single trials');
end

% to accommodate old behavior, only use trials that are completely filled
% with data
if cfg.vartrllength==1
  usetrial = false(size(data.trial),1);
  for k = 1:size(data.trial,1)
    usetrial(k,1) = all(isfinite(reshape(data.trial(k,:,:),[],1)));
  end
  fprintf('removing %d trials that contain NaNs\n', size(data.trial,1)-sum(usetrial));
  tmpcfg        = [];
  tmpcfg.trials = find(usetrial);
  data          = ft_selectdata(tmpcfg, data);% restore the provenance information
  [cfg, data]   = rollback_provenance(cfg, data);
end

[nrpt, nchan, nsmp] = size(data.trial);
if ~keeptrials
  avg = reshape(nanmean(data.trial,1),       [nchan nsmp]);
  dof = reshape(sum(isfinite(data.trial),1), [nchan nsmp]);
  var = reshape(nanvar(data.trial,0,1),      [nchan nsmp]);
else
  % nothing required here
end  

% compute the covariance matrix, if requested
if computecov
  if ischar(cfg.covariancewindow)
    switch cfg.covariancewindow
      case 'prestim'
        cfg.covariancewindow = [-inf 0];
      case 'poststim'
        cfg.covariancewindow = [0 inf];
      case 'all'
        cfg.covariancewindow = data.time([1 end]);
      otherwise
        ft_error('unsupported specification of cfg.covariancewindow');
    end
  end
  
  tmpcfg         = [];
  tmpcfg.latency = cfg.covariancewindow;
  tmpdata        = ft_selectdata(tmpcfg, data);
  
  if cfg.vartrllength==1
    % a second round of trial selection is needed here to accommodate old
    % behavior
    usetrial = false(nrpt,1);
    for k = 1:nrpt
      usetrial(k,1) = all(isfinite(reshape(tmpdata.trial(k,:,:),[],1)));
    end
    fprintf('removing %d trials that contain NaNs for the covariance computation\n', nrpt-sum(usetrial));
    tmpcfg        = [];
    tmpcfg.trials = find(usetrial);
    tmpdata       = ft_selectdata(tmpcfg, tmpdata);
  end
  
  % pre-allocate memory space for the covariance matrices
  if keeptrials
    covsig = nan(nrpt, nchan, nchan);
  else
    covsig = zeros(nchan, nchan);
    allsmp = 0;
  end
  
  % compute the covariance per trial
  for k = 1:nrpt
    dat    = reshape(tmpdata.trial(k,:,:), [nchan nsmp]);
    datsmp = isfinite(dat);
    if ~all(ismember(sum(datsmp), [0 nchan]))
      ft_error('channel specific NaNs are not supported for covariance computation');
    end
    numsmp = sum(datsmp(1,:));
    if istrue(cfg.removemean)
      dat  = ft_preproc_baselinecorrect(dat);
      smp = max(numsmp-1,1);
    end
    dat(~datsmp)  = 0;

    if keeptrials
      covsig(k,:,:) = dat*dat'./numsmp;
    else
      covsig = covsig + dat*dat';
      allsmp = allsmp + numsmp;
      % normalisation will be done after the for-loop
    end
  end
  
  if ~keeptrials 
    covsig = covsig./allsmp;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timelock = keepfields(data, {'time' 'grad', 'elec', 'opto', 'topo', 'topolabel', 'unmixing'});
if ~keeptrials
  timelock.avg        = avg;
  timelock.var        = var;
  timelock.dof        = dof;
  timelock.label      = data.label(:);
  timelock.dimord     = 'chan_time';
else
  timelock        = copyfields(data, timelock, {'trial' 'sampleinfo', 'trialinfo'});
  timelock.dimord = 'rpt_chan_time';
end
if computecov
  timelock.cov = covsig;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance timelock
ft_postamble history    timelock
ft_postamble savevar    timelock
