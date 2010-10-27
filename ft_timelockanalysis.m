function [timelock] = ft_timelockanalysis(cfg, data)

% FT_TIMELOCKANALYSIS performs timelocked analysis such as averaging
% and covariance computation
%
% [timelock] = ft_timelockanalysis(cfg, data)
%
% The data should be organised in a structure as obtained from the
% FT_PREPROCESSING function. The configuration should be according to
%
%   cfg.channel            = Nx1 cell-array with selection of channels (default = 'all'),
%                            see FT_CHANNELSELECTION for details
%   cfg.trials             = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.latency            = [begin end] in seconds, or 'minperlength', 'maxperlength', 
%                            'prestim', 'poststim' (default = 'maxperlength')
%   cfg.covariance         = 'no' or 'yes'
%   cfg.covariancewindow   = [begin end]
%   cfg.blcovariance       = 'no' or 'yes'
%   cfg.blcovariancewindow = [begin end]
%   cfg.keeptrials         = 'yes' or 'no', return individual trials or average (default = 'no')
%   cfg.normalizevar       = 'N' or 'N-1' (default = 'N-1')
%   cfg.removemean         = 'no' or 'yes' for covariance computation (default = 'yes')
%   cfg.vartrllength       = 0, 1 or 2 (see below)
%
% Depending on cfg.vartrllength, variable trials and missing values
% are treated differently:
%   0 - do not accept variable length trials [default]
%   1 - accept variable length trials, but only take those trials in which
%       data is present in both the average and the covariance window
%   2 - accept variable length trials, use all available trials
%       the available samples in every trial will be used for the
%       average and covariance computation. Missing values are replaced
%       by NaN and are not included in the computation.
%
% See also FT_TIMELOCKGRANDAVERAGE, FT_TIMELOCKSTATISTICS

% FIXME if input is one raw trial, the covariance is not computed correctly
% 
% Undocumented local options:
% cfg.feedback
% cfg.normalizecov
% cfg.preproc
% cfg.inputfile  = one can specifiy preanalysed saved data as input
% cfg.outputfile = one can specify output as file to save to disk

% This function depends on PREPROC which has the following options:
% cfg.absdiff
% cfg.blc
% cfg.blcwindow
% cfg.boxcar
% cfg.bpfilter
% cfg.bpfiltord
% cfg.bpfilttype
% cfg.bpfreq
% cfg.derivative
% cfg.detrend
% cfg.dftfilter
% cfg.dftfreq
% cfg.hilbert
% cfg.hpfilter
% cfg.hpfiltord
% cfg.hpfilttype
% cfg.hpfreq
% cfg.implicitref
% cfg.lpfilter
% cfg.lpfiltord
% cfg.lpfilttype
% cfg.lpfreq
% cfg.medianfilter
% cfg.medianfiltord
% cfg.rectify
% cfg.refchannel
% cfg.reref

% Copyright (C) 2003-2006, Markus Bauer
% Copyright (C) 2003-2006, Robert Oostenveld
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

fieldtripdefs

% set the defaults
if ~isfield(cfg, 'channel'),            cfg.channel = 'all';                     end
if ~isfield(cfg, 'trials'),             cfg.trials = 'all';                      end
if ~isfield(cfg, 'keeptrials'),         cfg.keeptrials = 'no';                   end
if ~isfield(cfg, 'latency'),            cfg.latency = 'maxperlength';            end
if ~isfield(cfg, 'covariance'),         cfg.covariance = 'no';                   end
if ~isfield(cfg, 'blcovariance'),       cfg.blcovariance = 'no';                 end
if ~isfield(cfg, 'removemean'),         cfg.removemean = 'yes';                  end
if ~isfield(cfg, 'vartrllength'),       cfg.vartrllength = 0;                    end
if ~isfield(cfg, 'feedback'),           cfg.feedback = 'text';                   end
if ~isfield(cfg, 'inputfile'),          cfg.inputfile = [];                      end
if ~isfield(cfg, 'outputfile'),         cfg.outputfile = [];                     end

hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    data = loadvar(cfg.inputfile, 'data');
  end
end

% check if the input data is valid for this function
data = checkdata(data, 'datatype', {'raw', 'comp'}, 'feedback', 'yes', 'hastrialdef', 'yes', 'hasoffset', 'yes');

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'deprecated',  {'normalizecov', 'normalizevar'});

% convert average to raw data for convenience, the output will be an average again
% the purpose of this is to allow for repeated baseline correction, filtering and other preproc options that timelockanalysis supports
data = data2raw(data);

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data = ft_selectdata(data, 'rpt', cfg.trials);  
end

ntrial = length(data.trial);

% ensure that the preproc specific options are located in the cfg.preproc substructure
cfg = checkconfig(cfg, 'createsubcfg',  {'preproc'});

% preprocess the data, i.e. apply filtering, baselinecorrection, etc.
fprintf('applying preprocessing options\n');
for i=1:ntrial
  [data.trial{i}, data.label, data.time{i}, cfg.preproc] = preproc(data.trial{i}, data.label, data.fsample, cfg.preproc, data.offset(i));
end

% determine the channels of interest
cfg.channel = ft_channelselection(cfg.channel, data.label);
chansel     = match_str(data.label, cfg.channel);
nchan       = length(cfg.channel);  % number of channels
numsamples  = zeros(ntrial,1);      % number of selected samples in each trial, is determined later

% determine the duration of each trial
for i=1:ntrial
  begsamplatency(i) = min(data.time{i});
  endsamplatency(i) = max(data.time{i});
end

% automatically determine the latency window which is possible in all trials
minperlength = [max(begsamplatency) min(endsamplatency)];
maxperlength = [min(begsamplatency) max(endsamplatency)];
maxtrllength = round((max(endsamplatency)-min(begsamplatency))*data.fsample) + 1;       % in samples
abstimvec    = ([1:maxtrllength] + min(data.offset) -1)./data.fsample;                  % in seconds

% latency window for averaging and variance computation is given in seconds
if (strcmp(cfg.latency, 'minperlength'))
  cfg.latency = [];
  cfg.latency(1) = minperlength(1);
  cfg.latency(2) = minperlength(2);
end
if (strcmp(cfg.latency, 'maxperlength'))
  cfg.latency = [];
  cfg.latency(1) = maxperlength(1);
  cfg.latency(2) = maxperlength(2);
end
if (strcmp(cfg.latency, 'prestim'))
  cfg.latency = [];
  cfg.latency(1) = maxperlength(1);
  cfg.latency(2) = 0;
end
if (strcmp(cfg.latency, 'poststim'))
  cfg.latency = [];
  cfg.latency(1) = 0;
  cfg.latency(2) = maxperlength(2);
end

% check whether trial type (varlength) matches the user-specified type
switch cfg.vartrllength
  case 0
    if ~all(minperlength==maxperlength)
      error('data has variable trial lengths, you specified not to accept that !');
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
    error('unknown value for vartrllength');
end

% check whether the time window fits with the data
if (cfg.latency(1) < maxperlength(1))  cfg.latency(1) = maxperlength(1);
  warning('Correcting begin latency of averaging window');
end
if (cfg.latency(2) > maxperlength(2))  cfg.latency(2) = maxperlength(2);
  warning('Correcting end latency of averaging window');
end
if cfg.latency(1)>cfg.latency(2)
  error('invalid latency window specified');
end

if strcmp(cfg.covariance, 'yes')
  if ~isfield(cfg, 'covariancewindow')
    % this used to be by default 'poststim', but that is not ideal as default
    error('the option cfg.covariancewindow is required');
  end
  % covariance window is given in seconds
  if (strcmp(cfg.covariancewindow, 'minperlength'))
    cfg.covariancewindow = [];
    cfg.covariancewindow(1) = minperlength(1);
    cfg.covariancewindow(2) = minperlength(2);
  end
  if (strcmp(cfg.covariancewindow, 'maxperlength'))
    cfg.covariancewindow = [];
    cfg.covariancewindow(1) = maxperlength(1);
    cfg.covariancewindow(2) = maxperlength(2);
  end
  if (strcmp(cfg.covariancewindow, 'prestim'))
    cfg.covariancewindow = [];
    cfg.covariancewindow(1) = maxperlength(1);
    cfg.covariancewindow(2) = 0;
  end
  if (strcmp(cfg.covariancewindow, 'poststim'))
    cfg.covariancewindow = [];
    cfg.covariancewindow(1) = 0;
    cfg.covariancewindow(2) = maxperlength(2);
  end
  % check whether the time window fits with the data
  if (cfg.covariancewindow(1) < maxperlength(1))
    cfg.covariancewindow(1) = maxperlength(1);
    warning('Correcting begin latency of covariance window');
  end
  if (cfg.covariancewindow(2) > maxperlength(2))
    cfg.covariancewindow(2) = maxperlength(2);
    warning('Correcting end latency of covariance window');
  end
  if cfg.covariancewindow(1)==cfg.covariancewindow(2)
    error('Cannot compute covariance over a window of only one sample');
  end
  if cfg.covariancewindow(1)>cfg.covariancewindow(2)
    error('Cannot compute covariance over negative timewindow');
  end
end

if strcmp(cfg.blcovariance, 'yes')
  if ~isfield(cfg, 'blcovariancewindow')
    % this used to be by default 'prestim', but that is not ideal as default
    error('the option cfg.blcovariancewindow is required');
  end
  % covariance window is given in seconds
  if (strcmp(cfg.blcovariancewindow, 'minperlength'))
    cfg.blcovariancewindow = [];
    cfg.blcovariancewindow(1) = minperlength(1);
    cfg.blcovariancewindow(2) = minperlength(2);
  end
  if (strcmp(cfg.blcovariancewindow, 'maxperlength'))
    cfg.blcovariancewindow = [];
    cfg.blcovariancewindow(1) = maxperlength(1);
    cfg.blcovariancewindow(2) = maxperlength(2);
  end
  if (strcmp(cfg.blcovariancewindow, 'prestim'))
    cfg.blcovariancewindow = [];
    cfg.blcovariancewindow(1) = maxperlength(1);
    cfg.blcovariancewindow(2) = 0;
  end
  if (strcmp(cfg.blcovariancewindow, 'poststim'))
    cfg.blcovariancewindow = [];
    cfg.blcovariancewindow(1) = 0;
    cfg.blcovariancewindow(2) = maxperlength(2);
  end
  % check whether the time window fits with the data
  if (cfg.blcovariancewindow(1) < maxperlength(1))
    cfg.blcovariancewindow(1) = maxperlength(1);
    warning('Correcting begin latency of covariance window');
  end
  if (cfg.blcovariancewindow(2) > maxperlength(2))
    cfg.blcovariancewindow(2) = maxperlength(2);
    warning('Correcting end latency of covariance window');
  end
  if cfg.blcovariancewindow(1)==cfg.blcovariancewindow(2)
    error('Cannot compute covariance over a window of only one sample');
  end
  if cfg.blcovariancewindow(1)>cfg.blcovariancewindow(2)
    error('Cannot compute covariance over negative timewindow');
  end
end

% pre-allocate some memory space for the covariance matrices
if strcmp(cfg.covariance, 'yes')
  covsig = nan*zeros(ntrial, nchan, nchan);
  numcovsigsamples = zeros(ntrial,1);
end
if strcmp(cfg.blcovariance, 'yes')
  covbl = nan*zeros(ntrial, nchan, nchan);
  numcovblsamples = zeros(ntrial,1);
end

begsampl = nearest(abstimvec, cfg.latency(1));
endsampl = nearest(abstimvec, cfg.latency(2));
maxwin   = endsampl-begsampl+1;
s        = zeros(nchan, maxwin);    % this will contain the sum
ss       = zeros(nchan, maxwin);    % this will contain the squared sum
dof      = zeros(1, maxwin);
if (strcmp(cfg.keeptrials,'yes'))
  singtrial = nan*zeros(ntrial, nchan, maxwin);
end

progress('init', cfg.feedback, 'averaging trials');
% do all the computations
for i=1:ntrial
  % fprintf('averaging trial %d of %d\n', i, ntrial);
  progress(i/ntrial, 'averaging trial %d of %d\n', i, ntrial);

  % determine whether the data in this trial can be used for all the requested computations
  switch cfg.vartrllength
    case 0
      % include this trial in any case since validation of data already done
      usetrial = 1;
    case 1
      % include this trial only if the data are complete in all specified windows
      usetrial = 1;
      if (begsamplatency(i)>cfg.latency(1) || endsamplatency(i)<cfg.latency(2))
        usetrial = 0;
      elseif strcmp(cfg.covariance,'yes') && (begsamplatency(i)>cfg.covariancewindow(1) || endsamplatency(i)<cfg.covariancewindow(2))
        usetrial = 0;
      elseif strcmp(cfg.blcovariance,'yes') && (begsamplatency(i)>cfg.blcovariancewindow(1) || endsamplatency(i)<cfg.blcovariancewindow(2))
        usetrial = 0;
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
  if (begsamplatency(i) <= cfg.latency(2)) && (endsamplatency(i) >= cfg.latency(1))
    begsampl = nearest(data.time{i}, cfg.latency(1));
    endsampl = nearest(data.time{i}, cfg.latency(2));
    numsamples(i) = endsampl-begsampl+1;
    if (cfg.latency(1)<begsamplatency(i))
      trlshift =round((begsamplatency(i)-cfg.latency(1))*data.fsample);
    else
      trlshift = 0;
    end
    windowsel = (1:numsamples(i))+trlshift;
    dat = data.trial{i}(chansel, begsampl:endsampl);
    if (strcmp(cfg.keeptrials,'yes'))
      % do not add the padded zeros to the 3D array, but keep the NaNs there to indicate missing values
      singtrial(i,:,:) = [nan*zeros(nchan, trlshift) dat nan*zeros(nchan,(maxwin-numsamples(i)-trlshift))];
    end
    dat = [zeros(nchan, trlshift) dat zeros(nchan,(maxwin-numsamples(i)-trlshift))];
    s  = s  + dat;            % compute the sum
    ss = ss + dat.^2;         % compute the squared sum
    % count the number of samples that went into the sum
    dof(windowsel) = dof(windowsel) + 1;
    usetrial = 1; % to indicate that this trial could be used
  end

  if strcmp(cfg.covariance, 'yes')
    begsampl = nearest(data.time{i}, cfg.covariancewindow(1));
    endsampl = nearest(data.time{i}, cfg.covariancewindow(2));
    numcovsigsamples(i) = endsampl-begsampl+1;
    % select the relevant samples from this trial, do NOT pad with zeros
    dat = data.trial{i}(chansel, begsampl:endsampl);
    if ~isempty(dat)  % we did not exlude this case above
      if strcmp(cfg.removemean, 'yes')
        dat = ft_preproc_baselinecorrect(dat);
      end
      covsig(i,:,:) = dat * dat';
    end
  end

  if strcmp(cfg.blcovariance, 'yes')
    begsampl = nearest(data.time{i}, cfg.blcovariancewindow(1));
    endsampl = nearest(data.time{i}, cfg.blcovariancewindow(2));
    numcovblsamples(i) = endsampl-begsampl+1;
    % select the relevant samples from this trial, do NOT pad with zeros
    dat = data.trial{i}(chansel, begsampl:endsampl);
    if ~isempty(dat)  % we did not exlude this case above
      if strcmp(cfg.removemean, 'yes')
        dat = ft_preproc_baselinecorrect(dat);
      end
      covbl(i,:,:) = dat * dat';
    end
  end
end % for ntrial
progress('close');

% compute the average
if ~any(numsamples)
  warning('no samples found in the specified time window, check option for vartrllength');
end
avg = s ./ repmat(dof(:)', [nchan 1]);


% compute the variance
% tmp1 = repmat(dof(:)', [nchan 1]);
% tmp2 = repmat(dof(:)', [nchan 1])-1;
% var = (ss - (s.^2)./tmp1) ./ tmp2;
dof = repmat(dof(:)', [nchan 1]);
var = (ss - (s.^2)./dof) ./ (dof-1);

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
      covsig = squeeze(nan_sum(covsig, 1)) / (sum(numcovsigsamples)-ntrial);
    else
      covsig = squeeze(nan_sum(covsig, 1)) / sum(numcovsigsamples);
    end
  end
end

% normalize the baseline covariance over all trials by the total number of samples in all trials
if strcmp(cfg.blcovariance, 'yes')
  if strcmp(cfg.keeptrials,'yes')
    for i=1:ntrial
      if strcmp(cfg.removemean, 'yes')
        covbl(i,:,:) = covbl(i,:,:) / (numcovblsamples(i)-1);
      else
        covbl(i,:,:) = covbl(i,:,:) / numcovblsamples(i);
      end
    end
  else
    if strcmp(cfg.removemean, 'yes')
      covbl = squeeze(sum(covbl, 1)) / (sum(numcovblsamples)-ntrial);
    else
      covbl = squeeze(sum(covbl, 1)) / sum(numcovblsamples);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timelock.avg        = avg;
timelock.var        = var;
timelock.fsample    = data.fsample;
timelock.time       = linspace(cfg.latency(1), cfg.latency(2), size(avg,2));
timelock.dof        = dof;
timelock.label      = data.label(chansel);
if (strcmp(cfg.keeptrials,'yes'))
  timelock.trial = singtrial;
  timelock.dimord = 'rpt_chan_time';
else
  timelock.dimord = 'chan_time';
end
if strcmp(cfg.covariance, 'yes')
  timelock.cov = covsig;
end
if strcmp(cfg.blcovariance, 'yes')
  timelock.blcov = covbl;
end
if isfield(data, 'grad')
  % copy the gradiometer array along
  timelock.grad = data.grad;
end
if isfield(data, 'elec')
  % copy the electrode array along
  timelock.elec = data.elec;
end
if isfield(data, 'trialinfo') && strcmp(cfg.keeptrials, 'yes')
  % copy the trialinfo into the output
  % but not the sampleinfo
  timelock.trialinfo = data.trialinfo;
end

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id$';

% remember the configuration details of the input data
try cfg.previous = data.cfg; end
% remember the exact configuration details in the output 
timelock.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', timelock); % use the variable name "data" in the output file
end
