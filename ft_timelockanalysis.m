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
%   cfg.latency            = [begin end] in seconds, or 'all', 'minperiod', 'maxperiod',
%                            'prestim', 'poststim' (default = 'all')
%   cfg.covariance         = 'no' or 'yes' (default = 'no')
%   cfg.covariancewindow   = [begin end] in seconds, or 'all', 'minperiod', 'maxperiod',
%                            'prestim', 'poststim' (default = 'all')
%   cfg.keeptrials         = 'yes' or 'no', return individual trials or average (default = 'no')
%   cfg.removemean         = 'no' or 'yes' for covariance computation (default = 'yes')
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
% cfg.blcovariance
% cfg.blcovariancewindow
% cfg.normalizevar
% cfg.normalizecov
% cfg.vartrllength

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
cfg = ft_checkconfig(cfg, 'forbidden',  {'blcovariance', 'blcovariancewindow'});
cfg = ft_checkconfig(cfg, 'renamed',    {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamed',    {'blcwindow', 'baselinewindow'});

% set the defaults
cfg.latency      = ft_getopt(cfg, 'latency',     'all');
cfg.trials       = ft_getopt(cfg, 'trials',      'all', 1);
cfg.channel      = ft_getopt(cfg, 'channel',     'all');
cfg.keeptrials   = ft_getopt(cfg, 'keeptrials',  'no');
cfg.vartrllength = ft_getopt(cfg, 'vartrllength', 0);
cfg.feedback     = ft_getopt(cfg, 'feedback',     'text');
cfg.preproc      = ft_getopt(cfg, 'preproc',      []);
cfg.covariance       = ft_getopt(cfg, 'covariance',      'no');
cfg.covariancewindow = ft_getopt(cfg, 'covariancewindow', 'all');
cfg.removemean       = ft_getopt(cfg, 'removemean',      'yes');

% create logical flags for convenience
keeptrials = istrue(cfg.keeptrials);
computecov = istrue(cfg.covariance);

% ensure that the preproc specific options are located in the cfg.preproc substructure
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});
if ~isempty(cfg.preproc)
  % preprocess the data, i.e. apply filtering, baselinecorrection, etc.
  fprintf('applying preprocessing options\n');
  if ~isfield(cfg.preproc, 'feedback')
    cfg.preproc.feedback = cfg.feedback;
  end
  data = ft_preprocessing(cfg.preproc, data);
  [cfg.preproc, data] = rollback_provenance(cfg.preproc, data);
end

% compute the covariance matrix, if requested
if computecov
  tmpcfg = keepfields(cfg, {'trials', 'channel', 'tolerance', 'showcallinfo'});
  tmpcfg.latency = cfg.covariancewindow;
  datacov = ft_selectdata(tmpcfg, data);
  % restore the provenance information
  [~, datacov] = rollback_provenance(cfg, datacov); % not sure what to do here
  datacov      = ft_checkdata(datacov, 'datatype', 'timelock');
  
  if isfield(datacov, 'trial')
    [nrpt, nchan, nsmp] = size(datacov.trial);
  else
    % if the data structure has only a single trial
    nrpt = 1;
    [nchan, nsmp] = size(datacov.avg);
    datacov.trial = shiftdim(datacov.avg, -1);
    datacov       = rmfield(datacov, 'avg');
    datacov.dimord = 'rpt_chan_time';
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
    dat    = reshape(datacov.trial(k,:,:), [nchan nsmp]);
    datsmp = isfinite(dat);
    if ~all(ismember(sum(datsmp,1), [0 nchan]))
      ft_error('channel specific NaNs are not supported for covariance computation');
    end
    numsmp = sum(datsmp(1,:));
    if istrue(cfg.removemean)
      dat  = ft_preproc_baselinecorrect(dat);
      numsmp = max(numsmp-1,1);
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

% select trials and channels of interest
tmpcfg = keepfields(cfg, {'trials', 'channel', 'tolerance', 'latency', 'showcallinfo'});
data   = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

% convert to a timelock structure with trials kept and NaNs for missing
% data points, when there's only a single trial in the input data
% structure, this leads to an 'avg' field, rather than a 'trial' field,
% and also the trialinfo is removed, so keep separate before conversion
if isfield(data, 'trialinfo'), trialinfo = data.trialinfo; end
data = ft_checkdata(data, 'datatype', {'timelock+comp' 'timelock'});

if ~keeptrials && isfield(data, 'trial')
  [nrpt, nchan, nsmp] = size(data.trial);
  avg = reshape(nanmean(data.trial,1),       [nchan nsmp]);
  dof = reshape(sum(isfinite(data.trial),1), [nchan nsmp]);
  var = reshape(nanvar(data.trial,0,1),      [nchan nsmp]);
elseif ~keeptrials && ~isfield(data, 'trial')
  avg = data.avg;
  var = nan(size(data.avg));
  dof = double(isfinite(data.avg));
elseif keeptrials && isfield(data, 'trial')
  % nothing required here
elseif keeptrials && ~isfield(data, 'trial')
  % don't know whether this is a use case
  data.trial = shiftdim(data.avg, -1);
  if exist('trialinfo', 'var')
    data.trialinfo = trialinfo;
  end
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timelock = keepfields(data, {'time' 'grad', 'elec', 'opto', 'topo', 'topodimord', 'topolabel', 'unmixing', 'unmixingdimord', 'label'});
if ~keeptrials
  timelock.avg        = avg;
  timelock.var        = var;
  timelock.dof        = dof;
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
