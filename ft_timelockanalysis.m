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
%   cfg.covariance         = 'no' or 'yes'
%   cfg.covariancewindow   = [begin end]
%   cfg.keeptrials         = 'yes' or 'no', return individual trials or average (default = 'no')
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
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
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

% This function depends on PREPROC which has the following options:
% cfg.absdiff
% cfg.boxcar
% cfg.bpfilter
% cfg.bpfiltord
% cfg.bpfilttype
% cfg.bpfreq
% cfg.demean
% cfg.baselinewindow
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble distribute
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar data

% return immediately after distributed execution
if ~isempty(ft_getopt(cfg, 'distribute'))
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'raw', 'comp'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'deprecated',  {'normalizecov', 'normalizevar'});
cfg = ft_checkconfig(cfg, 'deprecated',  {'latency', 'blcovariance', 'blcovariancewindow'});
cfg = ft_checkconfig(cfg, 'renamed',     {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamed',     {'blcwindow', 'baselinewindow'});

% set the defaults
if ~isfield(cfg, 'channel'),       cfg.channel      = 'all';  end
if ~isfield(cfg, 'trials'),        cfg.trials       = 'all';  end
if ~isfield(cfg, 'keeptrials'),    cfg.keeptrials   = 'no';   end
if ~isfield(cfg, 'covariance'),    cfg.covariance   = 'no';   end
if ~isfield(cfg, 'removemean'),    cfg.removemean   = 'yes';  end
if ~isfield(cfg, 'vartrllength'),  cfg.vartrllength = 0;      end
if ~isfield(cfg, 'feedback'),      cfg.feedback     = 'text'; end
if ~isfield(cfg, 'preproc'),       cfg.preproc      = [];     end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data = ft_selectdata(data, 'rpt', cfg.trials);  
end

ntrial = length(data.trial);

% ensure that the preproc specific options are located in the cfg.preproc substructure
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});

if ~isempty(cfg.preproc)
  % preprocess the data, i.e. apply filtering, baselinecorrection, etc.
  fprintf('applying preprocessing options\n');
  data = ft_preprocessing(cfg.preproc, data);
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
  offset(i)         = time2offset(data.time{i}, data.fsample);
end

% automatically determine the latency window which is possible in all trials
minperlength = [max(begsamplatency) min(endsamplatency)];
maxperlength = [min(begsamplatency) max(endsamplatency)];
maxtrllength = round((max(endsamplatency)-min(begsamplatency))*data.fsample) + 1;       % in samples
abstimvec    = ([1:maxtrllength] + min(offset) -1)./data.fsample;                  % in seconds

latency      = [];
latency(1)   = maxperlength(1);
latency(2)   = maxperlength(2);

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

if strcmp(cfg.covariance, 'yes')
  if ~isfield(cfg, 'covariancewindow')
    warning('the option cfg.covariancewindow is not specified, taking all time points');
    cfg.covariancewindow = latency;
  end
  if ischar(cfg.covariancewindow)
    switch cfg.covariancewindow
    case 'prestim'
      cfg.covariancewindow = [latency(1) 0];
    case 'poststim'
      cfg.covariancewindow = [0 latency(2)];
    case 'minperlength'
      error('cfg.covariancewindow = ''minperlength'' is not supported anymore');
    case 'maxperlength'
      error('cfg.covariancewindow = ''maxperlength'' is not supported anymore');
    otherwise
      error('unsupported specification of cfg.covariancewindow');
    end
  end
end

% pre-allocate some memory space for the covariance matrices
if strcmp(cfg.covariance, 'yes')
  covsig = nan*zeros(ntrial, nchan, nchan);
  numcovsigsamples = zeros(ntrial,1);
end

begsampl = nearest(abstimvec, latency(1));
endsampl = nearest(abstimvec, latency(2));
maxwin   = endsampl-begsampl+1;
s        = zeros(nchan, maxwin);    % this will contain the sum
ss       = zeros(nchan, maxwin);    % this will contain the squared sum
dof      = zeros(1, maxwin);
if (strcmp(cfg.keeptrials,'yes'))
  singtrial = nan*zeros(ntrial, nchan, maxwin);
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
        warning(['trial ' num2str(i) ' not used for avg computation because it was not used for covariance computation']);
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
    if (latency(1)<begsamplatency(i))
      trlshift =round((begsamplatency(i)-latency(1))*data.fsample);
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

end % for ntrial
ft_progress('close');

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

if (dof > 1)
  var = (ss - (s.^2)./dof) ./ (dof-1);
else
  var = zeros(size(avg));
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
      covsig = squeeze(nansum(covsig, 1)) / (sum(numcovsigsamples)-ntrial);
    else
      covsig = squeeze(nansum(covsig, 1)) / sum(numcovsigsamples);
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

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous data
ft_postamble history timelock
ft_postamble savevar timelock

