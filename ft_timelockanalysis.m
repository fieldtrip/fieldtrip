function [timelock] = ft_timelockanalysis(cfg, data)

% FT_TIMELOCKANALYSIS computes the timelocked average ERP/ERF and
% computes the covariance matrix
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
cfg.covariance   = ft_getopt(cfg, 'covariance',  'no');
cfg.removemean   = ft_getopt(cfg, 'removemean',  'yes');
cfg.vartrllength = ft_getopt(cfg, 'vartrllength', 0);
cfg.feedback     = ft_getopt(cfg, 'feedback',     'text');
cfg.preproc      = ft_getopt(cfg, 'preproc',      []);

% ensure that the preproc specific options are located in the cfg.preproc substructure
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});

% select trials of interest
tmpcfg = [];
tmpcfg.trials = cfg.trials;
tmpcfg.channel = cfg.channel;
data = ft_selectdata(tmpcfg, data);
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

% determine the size of the data
ntrial      = length(data.trial);
nchan       = length(data.label);   % number of channels
numsamples  = zeros(ntrial,1);      % number of selected samples in each trial, is determined later

if ntrial==0, ft_error('Number of trials selected in data is zero');   end
if nchan==0,  ft_error('Number of channels selected in data is zero'); end

% determine the duration of each trial
begsamplatency = zeros(1,ntrial);
endsamplatency = zeros(1,ntrial);
offset         = zeros(1,ntrial);
for i=1:ntrial
  begsamplatency(i) = min(data.time{i});
  endsamplatency(i) = max(data.time{i});
  offset(i)         = time2offset(data.time{i}, data.fsample);
end

% automatically determine the latency window which is possible in all trials
minperlength = [max(begsamplatency) min(endsamplatency)];
maxperlength = [min(begsamplatency) max(endsamplatency)];
maxtrllength = round((max(endsamplatency)-min(begsamplatency))*data.fsample) + 1;       % in samples
abstimvec    = ((1:maxtrllength) + min(offset) -1)./data.fsample;                       % in seconds

latency      = [];
latency(1)   = maxperlength(1);
latency(2)   = maxperlength(2);

% check whether trial type (varlength) matches the user-specified type
switch cfg.vartrllength
  case 0
    if ~all(minperlength==maxperlength)
      ft_error('data has variable trial lengths, you specified not to accept that');
    end
  case 1
    if all(minperlength==maxperlength)
      disp('data is of type fixed length !');
    end
  case 2
    if strcmp(cfg.keeptrials,'yes')
      disp('processing and keeping variable length single trials');
    end
  otherwise
    ft_error('unknown value for vartrllength');
end

if strcmp(cfg.covariance, 'yes')
  if ~isfield(cfg, 'covariancewindow')
    ft_warning('the option cfg.covariancewindow is not specified, taking all time points');
    cfg.covariancewindow = latency;
  end
  if ischar(cfg.covariancewindow)
    switch cfg.covariancewindow
      case 'prestim'
        cfg.covariancewindow = [latency(1) 0];
      case 'poststim'
        cfg.covariancewindow = [0 latency(2)];
      case 'all'
        cfg.covariancewindow = latency;
      case 'minperlength'
        ft_error('cfg.covariancewindow = ''minperlength'' is not supported anymore');
      case 'maxperlength'
        ft_error('cfg.covariancewindow = ''maxperlength'' is not supported anymore');
      otherwise
        ft_error('unsupported specification of cfg.covariancewindow');
    end
  end
end

% pre-allocate some memory space for the covariance matrices
if strcmp(cfg.covariance, 'yes')
  covsig = nan(ntrial, nchan, nchan);
  numcovsigsamples = zeros(ntrial,1);
end

begsampl = nearest(abstimvec, latency(1));
endsampl = nearest(abstimvec, latency(2));
maxwin   = endsampl-begsampl+1;
s        = zeros(nchan, maxwin);    % this will contain the sum
ss       = zeros(nchan, maxwin);    % this will contain the squared sum
dof      = zeros(1, maxwin);
if (strcmp(cfg.keeptrials,'yes'))
  singtrial = nan(ntrial, nchan, maxwin);
end

ft_progress('init', cfg.feedback, 'averaging trials');
% do all the computations
for i=1:ntrial
  % fprintf('averaging trial %d of %d\n', i, ntrial);
  ft_progress(i/ntrial, 'averaging trial %d of %d\n', i, ntrial);

  % determine whether the data in this trial can be used for all the requested computations
  switch cfg.vartrllength
    case 0
      % include this trial in any case since validation of data already done
      usetrial = 1;
    case 1
      % include this trial only if the data are complete in all specified windows
      usetrial = 1;
      %       if (begsamplatency(i)>latency(1) || endsamplatency(i)<latency(2))
      %         usetrial = 0;
      %       elseif strcmp(cfg.covariance,'yes') && (begsamplatency(i)>cfg.covariancewindow(1) || endsamplatency(i)<cfg.covariancewindow(2))
      if strcmp(cfg.covariance,'yes') && (begsamplatency(i)>cfg.covariancewindow(1) || endsamplatency(i)<cfg.covariancewindow(2))
        usetrial = 0;
        ft_warning(['trial ' num2str(i) ' not used for avg computation because it was not used for covariance computation']);
      end
    case 2
      % include this trial if any data points are present in any of the specified windows
      % this is handled automatically by the code below
      usetrial = 1;
  end

  if ~usetrial
    continue;
  end

  % for average and variance
  if (begsamplatency(i) <= latency(2)) && (endsamplatency(i) >= latency(1))
    begsampl = nearest(data.time{i}, latency(1));
    endsampl = nearest(data.time{i}, latency(2));
    numsamples(i) = endsampl-begsampl+1;
    dat = data.trial{i}(:, begsampl:endsampl);
    if isfield(data, 'sampleinfo')
      tmpsampl = data.sampleinfo(i,1):data.sampleinfo(i,2);
      data.sampleinfo(i,:) = tmpsampl([begsampl endsampl]);
    end
    
    if (latency(1)<begsamplatency(i))
      trlshift = floor((begsamplatency(i)-latency(1))*data.fsample);
    else
      trlshift = 0;
    end
    windowsel = (1:numsamples(i))+trlshift;
    if (strcmp(cfg.keeptrials,'yes'))
      % do not pad with zeros, but keep the NaNs to indicate missing values
      singtrial(i,:,windowsel) = dat;
    end
    s (:,windowsel) = s (:,windowsel) + dat;            % compute the sum
    ss(:,windowsel) = ss(:,windowsel) + dat.^2;         % compute the sum of squares
    % count the number of samples that went into the sum
    dof(windowsel) = dof(windowsel) + 1;
    usetrial = 1; % to indicate that this trial could be used
  end

  if strcmp(cfg.covariance, 'yes')
    begsampl = nearest(data.time{i}, cfg.covariancewindow(1));
    endsampl = nearest(data.time{i}, cfg.covariancewindow(2));
    numcovsigsamples(i) = endsampl-begsampl+1;
    % select the relevant samples from this trial, do NOT pad with zeros
    dat = data.trial{i}(:, begsampl:endsampl);
    if ~isempty(dat)  % we did not exlude this case above
      if strcmp(cfg.removemean, 'yes')
        dat = ft_preproc_baselinecorrect(dat);
      end
      covsig(i,:,:) = dat * dat';
    end
  end

end % for ntrial
ft_progress('close');

% compute the average
if ~any(numsamples)
  ft_warning('no samples found in the specified time window, check option for vartrllength');
end
avg = s ./ repmat(dof(:)', [nchan 1]);


% compute the variance
% tmp1 = repmat(dof(:)', [nchan 1]);
% tmp2 = repmat(dof(:)', [nchan 1])-1;
% var = (ss - (s.^2)./tmp1) ./ tmp2;
dof = repmat(dof(:)', [nchan 1]);

if any(dof(:) > 1)
  var = (ss - (s.^2)./dof) ./ (dof-1);
else
  var = nan(size(avg));
end

% normalize the covariance over all trials by the total number of samples in all trials
if strcmp(cfg.covariance, 'yes')
  if strcmp(cfg.keeptrials,'yes')
    for i=1:ntrial
      if strcmp(cfg.removemean, 'yes')
        covsig(i,:,:) = covsig(i,:,:) / (numcovsigsamples(i)-1);
      else
        covsig(i,:,:) = covsig(i,:,:) / numcovsigsamples(i);
      end
    end
  else
    if strcmp(cfg.removemean, 'yes')
      covsig = shiftdim(nansum(covsig, 1)) / (sum(numcovsigsamples)-ntrial);
    else
      covsig = shiftdim(nansum(covsig, 1)) / sum(numcovsigsamples);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timelock.avg        = avg;
timelock.var        = var;
%timelock.fsample    = data.fsample; % timelock.fsample is obsolete
timelock.time       = linspace(latency(1), latency(2), size(avg,2));
timelock.dof        = dof;
timelock.label      = data.label(:);
if (strcmp(cfg.keeptrials,'yes'))
  timelock.trial = singtrial;
  timelock.dimord = 'rpt_chan_time';
else
  timelock.dimord = 'chan_time';
end
if strcmp(cfg.covariance, 'yes')
  timelock.cov = covsig;
end

% some fields from the input should be copied over in the output
timelock = copyfields(data, timelock, {'grad', 'elec', 'opto', 'topo', 'topolabel', 'unmixing'});

if isfield(data, 'trialinfo') && strcmp(cfg.keeptrials, 'yes')
  timelock.trialinfo = data.trialinfo;
end
if isfield(data, 'sampleinfo') && strcmp(cfg.keeptrials, 'yes')
  timelock.sampleinfo = data.sampleinfo;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance timelock
ft_postamble history    timelock
ft_postamble savevar    timelock
