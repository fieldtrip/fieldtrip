function [timelock] = ft_timelockanalysis(cfg, data)

% FT_TIMELOCKANALYSIS computes the timelocked average ERP/ERF and optionally computes
% the covariance matrix over the specified time window.
%
% Use as
%   [timelock] = ft_timelockanalysis(cfg, data)
%
% The data should be organised in a structure as obtained from the FT_PREPROCESSING
% function. The configuration should be according to
%   cfg.channel            = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.latency            = [begin end] in seconds, or 'all', 'minperiod', 'maxperiod', 'prestim', 'poststim' (default = 'all')
%   cfg.trials             = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.keeptrials         = 'yes' or 'no', return individual trials or average (default = 'no')
%   cfg.nanmean            = string, can be 'yes' or 'no' (default = 'yes')
%   cfg.normalizevar       = 'N' or 'N-1' (default = 'N-1')
%   cfg.covariance         = 'no' or 'yes' (default = 'no')
%   cfg.covariancewindow   = [begin end] in seconds, or 'all', 'minperiod', 'maxperiod', 'prestim', 'poststim' (default = 'all')
%   cfg.removemean         = 'yes' or 'no', for the covariance computation (default = 'yes')
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

% FIXME if input is one raw trial, the covariance is not computed correctly
%
% Undocumented local options:
%   cfg.feedback
%   cfg.preproc
%
% Deprecated options:
%   cfg.blcovariance
%   cfg.blcovariancewindow
%   cfg.normalizecov
%   cfg.vartrllength

% Copyright (C) 2018, Jan-Mathijs Schoffelen
% Copyright (C) 2003-2006, Markus Bauer
% Copyright (C) 2003-2022, Robert Oostenveld
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

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels', 'trial'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'forbidden',  {'normalizecov'});
cfg = ft_checkconfig(cfg, 'forbidden',  {'blcovariance', 'blcovariancewindow'});
cfg = ft_checkconfig(cfg, 'renamed',    {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamed',    {'blcwindow', 'baselinewindow'});

% set the defaults
cfg.preproc           = ft_getopt(cfg, 'preproc'          , []);
cfg.channel           = ft_getopt(cfg, 'channel'          , 'all');
cfg.latency           = ft_getopt(cfg, 'latency'          , 'all');
cfg.trials            = ft_getopt(cfg, 'trials'           , 'all', 1);
cfg.keeptrials        = ft_getopt(cfg, 'keeptrials'       , 'no');
cfg.vartrllength      = ft_getopt(cfg, 'vartrllength'     , 0);
cfg.nanmean           = ft_getopt(cfg, 'nanmean'          , 'yes');
cfg.normalizevar      = ft_getopt(cfg, 'normalizevar'     , 'N-1');
cfg.covariance        = ft_getopt(cfg, 'covariance'       , 'no');
cfg.covariancewindow  = ft_getopt(cfg, 'covariancewindow' , 'all');
cfg.removemean        = ft_getopt(cfg, 'removemean'       , 'yes');
cfg.feedback          = ft_getopt(cfg, 'feedback'         , 'text');

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
  tmpcfg = keepfields(cfg, {'trials', 'channel', 'tolerance', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
  tmpcfg.latency = cfg.covariancewindow;
  datacov = ft_selectdata(tmpcfg, data);
  % restore the provenance information
  [dum, datacov] = rollback_provenance(cfg, datacov); % not sure what to do here
  datacov      = ft_checkdata(datacov, 'datatype', 'timelock');

  if isfield(datacov, 'trial')
    [nrpt, nchan, ntime] = size(datacov.trial);
  else
    % if the data structure has only a single trial
    nrpt = 1;
    [nchan, ntime] = size(datacov.avg);
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
    dat    = reshape(datacov.trial(k,:,:), [nchan ntime]);
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
orgcfg = cfg;
tmpcfg = keepfields(cfg, {'trials', 'channel', 'tolerance', 'latency', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
data   = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);
% do not use the default option returned by FT_SELECTDATA, but the original one for this function
cfg.nanmean = orgcfg.nanmean;

% do a sanity check
if isempty(data.trial)
  if ~isempty(cfg.trials)
    ft_error('there are no trials selected');
  else
    ft_error('there are no trials in the input data');
  end
end

if keeptrials
  % convert to a timelock structure with trials kept and NaNs for missing data points, when there's only a single trial in the input data
  % structure, this leads to an 'avg' field, rather than a 'trial' field, and also the trialinfo is removed, so keep separate before conversion
  if isfield(data, 'trialinfo'), trialinfo = data.trialinfo; end
  data = ft_checkdata(data, 'datatype', {'timelock+comp' 'timelock'});
  
  if keeptrials && isfield(data, 'trial')
    % nothing required here
  elseif keeptrials && ~isfield(data, 'trial')
    % don't know whether this is a use case
    data.trial = shiftdim(data.avg, -1);
    if exist('trialinfo', 'var')
      data.trialinfo = trialinfo;
    end
  end
  
elseif ~keeptrials
  % whether to normalize the variance with N or N-1, see VAR
  normalizewithN = strcmpi(cfg.normalizevar, 'N');
  
  % compute a running sum average/var etc. to save memory
  
  % the code below tries to construct a general time-axis where samples of all trials can fall on
  % find the earliest beginning and latest ending
  begtime = min(cellfun(@min, data.time));
  endtime = max(cellfun(@max, data.time));
  % find 'common' sampling rate
  fsample = 1./nanmean(cellfun(@mean, cellfun(@diff,data.time, 'uniformoutput', false)));
  % estimate number of samples
  nsmp = round((endtime-begtime)*fsample) + 1; % numerical round-off issues should be dealt with by this round, as they will/should never cause an extra sample to appear
  % construct general time-axis
  time = linspace(begtime, endtime, nsmp);
  
  nchan  = numel(data.label);
  ntrial = numel(data.trial);

  % placeholder for running sums
  tmpsum = zeros(nchan, length(time));
  tmpssq = tmpsum;
  tmpdof = tmpsum;
  
  begsmp = nan(ntrial, 1);
  endsmp = nan(ntrial, 1);

  % do a 2-pass running sum, sacrificing speed for numeric stability
  for i=1:ntrial
    begsmp(i) = nearest(time, data.time{i}(1));
    endsmp(i) = nearest(time, data.time{i}(end));
    
    tmp = data.trial{i};
    
    tmpdof(:,begsmp(i):endsmp(i)) = isfinite(tmp) + tmpdof(:,begsmp(i):endsmp(i));
    if istrue(cfg.nanmean)
      tmp(~isfinite(tmp)) = 0;
    end
    tmpsum(:,begsmp(i):endsmp(i)) = tmp    + tmpsum(:,begsmp(i):endsmp(i));
  end
  avgmat = tmpsum ./ tmpdof;

  tmpsum = zeros(nchan, length(time));
  for i=1:ntrial
    tmp = data.trial{i};
    
    tmp = tmp - avgmat(:,begsmp(i):endsmp(i));

    if istrue(cfg.nanmean)
      tmp(~isfinite(tmp)) = 0;
    end
    tmpsum(:,begsmp(i):endsmp(i)) = tmp    + tmpsum(:,begsmp(i):endsmp(i));
    tmpssq(:,begsmp(i):endsmp(i)) = tmp.^2 + tmpssq(:,begsmp(i):endsmp(i));
    
  end
  
  dofmat = tmpdof;
  %avgmat = tmpsum ./ tmpdof;
  varmat = tmpssq ./ tmpdof - (tmpsum ./ tmpdof).^2;
  if normalizewithN
    
    % just to be sure
    varmat(dofmat<=0) = NaN;
  else
    varmat = varmat .* (dofmat ./ (dofmat-1));
    
    % see https://stats.stackexchange.com/questions/4068/how-should-one-define-the-sample-variance-for-scalar-input
    % the fieldtrip/external/stats/nanvar implementation behaves differently here than Mathworks VAR and NANVAR implementations
    varmat(dofmat<=1) = NaN;
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timelock = keepfields(data, {'time', 'grad', 'elec', 'opto', 'topo', 'topodimord', 'topolabel', 'unmixing', 'unmixingdimord', 'label'});
if ~keeptrials
  timelock.avg        = avgmat;
  timelock.var        = varmat;
  timelock.dof        = dofmat;
  timelock.time       = time;
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
ft_postamble previous   data
ft_postamble provenance timelock
ft_postamble history    timelock
ft_postamble savevar    timelock
