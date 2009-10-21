function [timelock] = timelockanalysis(cfg, data)

% TIMELOCKANALYSIS performs timelocked analysis such as averaging
% and covariance computation
%
% [timelock] = timelockanalysis(cfg, data)
%
% The data should be organised in a structure as obtained from the
% PREPROCESSING function. The configuration should be according to
%
%   cfg.channel            = Nx1 cell-array with selection of channels (default = 'all'),
%                            see CHANNELSELECTION for details
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

% FIXME if input is one raw trial, the covariance is not computed correctly
% 
% Undocumented local options:
% cfg.feedback
% cfg.normalizecov
% cfg.preproc
%
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
% $Log: timelockanalysis.m,v $
% Revision 1.60  2009/03/23 21:21:16  roboos
% removed defaults for covariancewindow and blcovariancewindow, since the defaults were not optimal for most cases. Now the user is forced to specify the window.
%
% Revision 1.59  2009/02/04 16:44:06  roboos
% remove numsamples, numcovsamples and numblcovsamples from the output, since these are not used anywhere
%
% Revision 1.58  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.57  2009/01/12 13:05:20  sashae
% small change in call to checkconfig
%
% Revision 1.56  2008/11/11 18:59:26  sashae
% added call to checkconfig at end of function (trackconfig and checksize)
%
% Revision 1.55  2008/10/01 15:51:26  sashae
% call to checkconfig instead of createsubcfg
%
% Revision 1.54  2008/09/26 12:42:18  sashae
% checkconfig: checks if the input cfg is valid for this function
%
% Revision 1.53  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.52  2008/07/11 13:18:40  roboos
% removed all lnfilter references, added error to preprocessing and preproc
%
% Revision 1.51  2008/06/10 16:45:04  sashae
% replaced call to blc function with preproc_baselinecorrect
%
% Revision 1.50  2008/05/06 15:43:46  sashae
% change in trial selection, cfg.trials can be logical
%
% Revision 1.49  2008/01/29 18:06:13  sashae
% added option for trial selection
% removed some old code
% adding of offset field now done with checkdata (hasoffset='yes')
%
% Revision 1.48  2007/10/31 17:05:00  roboos
% represent dof as matrix instead of vector
%
% Revision 1.47  2007/05/30 11:41:10  roboos
% renamed dofvec into dof
%
% Revision 1.46  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.45  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.44  2006/10/04 07:10:08  roboos
% updated documentation
%
% Revision 1.43  2006/06/14 12:44:03  roboos
% removed the documentation for cfg.lnfilttype, since that option is not supported by preproc
%
% Revision 1.42  2006/06/14 11:53:54  roboos
% switched to using cfg.preproc substructure
%
% Revision 1.41  2006/06/13 14:48:11  ingnie
% updated documentation
%
% Revision 1.40  2006/06/06 16:57:51  ingnie
% updated documentation
%
% Revision 1.39  2006/04/25 17:06:28  ingnie
% updated documentation
%
% Revision 1.38  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.37  2006/02/01 12:26:00  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.36  2005/11/21 11:57:37  roboos
% replaced sum by nan_sum in computation of covariance for variable length data
%
% Revision 1.35  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.34  2005/06/17 10:56:12  roboos
% switched from local implementation of data conversion towards new data2raw private function
% this also fixed a problem that I observed with some simulated data, but I did not look into detail in that problem
%
% Revision 1.33  2005/05/17 17:50:39  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.32  2005/04/25 09:53:16  roboos
% fixed important bug in channel label ordering when specifying own cfg.channel
% fixed small bug for trials containing only one sample (these were left out)
% changed output cfg assignment of preprocessing cfg options
%
% Revision 1.31  2005/04/07 17:05:36  roboos
% replaced fprintf feedback by more flexible progress indicator
%
% Revision 1.30  2005/04/05 19:12:48  roboos
% added support for average input (converted to single trial) to allow for more interactive playing around with filter and blc settings for figures
% fixed error in preproc when cfg.implicitref was specified (implicitref was not included in the outpu of timelockanalysis)
%
% Revision 1.29  2005/02/16 08:44:11  roboos
% fixed weird typo in numsamples, don't know when that was introduced
%
% Revision 1.28  2004/12/09 17:33:06  roboos
% removed blc and bpfilter
% implemented filtering and all other preprocessing options with new private/preproc function
%
% Revision 1.27  2004/11/02 14:30:39  roboos
% changed indentation of the code, no functional changes
%
% Revision 1.26  2004/10/25 10:02:34  roboos
% fixed bug in covariance computation for 1 trial
% cleaned up the definition and detection of the length of the time axis in the input data
%
% Revision 1.25  2004/10/22 16:18:41  roboos
% also added dimord if keeptrials=no
%
% Revision 1.24  2004/09/22 10:20:27  roboos
% converted to use external subfunctions time2offset and offset2time
% and add offset field to data structure if it is missing
%
% Revision 1.23  2004/08/06 06:26:29  roboos
% only stylistic changes, nu functionality change
%
% Revision 1.22  2004/06/24 10:14:46  roberto
% implemented baselinecorrection (as undocumented/advanced feature)
% keep data.elec in output average just as data.grad
%
% Revision 1.21  2004/04/13 16:31:09  roberto
% fixed bug in dbstack selection of function filename for Matlab 6.1
%
% Revision 1.20  2004/04/13 14:25:24  roberto
% wrapped code-snippet around mfilename to make it compatible with Matlab 6.1
%
% Revision 1.19  2004/03/29 15:08:18  roberto
% added cfg.bpfilttype = fir|but to documentation
%
% Revision 1.18  2004/03/22 15:57:37  roberto
% restructured version output in configuration
%
% Revision 1.17  2004/03/22 15:47:23  roberto
% added version and previous cfg to output cfg
%
% Revision 1.16  2004/03/10 12:55:26  roberto
% fixed bug in output dimord, should be trial_chan_time instead of chan_time_trial
%
% Revision 1.15  2004/03/03 15:27:32  roberto
% incorporated couple of bug-fixes from Markus B. and one from Markus S.
% improved and restructured help
%
% Revision 1.14  2004/02/27 13:13:08  roberto
% fixed small error in the documentation
%
% Revision 1.13  2004/02/23 08:50:29  roberto
% added a dimord string to the output
%
% Revision 1.12  2004/02/17 15:36:47  roberto
% fixed multiple bugs that were still outstanding ue to the recent added support for variable trial length
%
% Revision 1.11  2004/02/17 09:30:05  roberto
% corrected automatic updating of latency windows
%
% Revision 1.10  2004/02/16 17:02:42  roberto
% added support for variable length trials, this required many changes
%
% Revision 1.9  2004/02/11 09:05:09  roberto
% removed default bandpass frequencies, changed default bpfiltord to
% 4 (from 25), added support for other filter types, removed the
% filtering options from the visible help
%
% Revision 1.8  2004/01/26 16:48:39  roberto
% added check for valid latency window
%
% Revision 1.7  2003/11/17 15:57:22  roberto
% fixed bug in counting number of channels in output data
%
% Revision 1.6  2003/11/17 15:15:02  roberto
% renamed default cfg item bpfltord into bpfiltord to fix code inconsistency
%
% Revision 1.5  2003/10/31 11:39:41  roberto
% switched from translate_channel_list to channelselection function
%
% Revision 1.4  2003/10/23 06:41:14  roberto
% added Markus' bandpass filtering
%
% Revision 1.3  2003/09/11 21:53:47  roberto
% switched to new channel group selection with translate_channel_list()
% fixed bug caused by removing too much code
%
% Revision 1.2  2003/04/23 10:14:32  roberto
% fixed some small bugs that I thought were fixed already earlier (!?)
%
% Revision 1.1.1.1  2003/04/17 12:35:19  roberto
% initial version under CVS control
%

fieldtripdefs

% check if the input data is valid for this function
data = checkdata(data, 'datatype', {'raw', 'comp'}, 'feedback', 'yes', 'hasoffset', 'yes');

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'deprecated',  {'normalizecov', 'normalizevar'});

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

% convert average to raw data for convenience, the output will be an average again
% the purpose of this is to allow for repeated baseline correction, filtering and other preproc options that timelockanalysis supports
data = data2raw(data);

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  if islogical(cfg.trials),  cfg.trials=find(cfg.trials);  end
  fprintf('selecting %d trials\n', length(cfg.trials));
  data.trial  = data.trial(cfg.trials);
  data.time   = data.time(cfg.trials);
  data.offset = data.offset(cfg.trials);
  % update the trial definition (trl)
  if isfield(data, 'cfg') % try to locate the trl in the nested configuration
    trl = findcfg(data.cfg, 'trl');
  else
    trl = [];
  end
  if isempty(trl)
    % a trial definition is expected in each continuous data set
    warning('could not locate the trial definition ''trl'' in the data structure');
  else
    cfg.trlold=trl;
    cfg.trl=trl(cfg.trials,:);
  end
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
cfg.channel = channelselection(cfg.channel, data.label);
chansel     = match_str(data.label, cfg.channel);
nchan       = length(cfg.channel);	% number of channels
numsamples  = zeros(ntrial,1);		% number of selected samples in each trial, is determined later

% determine the duration of each trial
for i=1:ntrial
  begsamplatency(i) = min(data.time{i});
  endsamplatency(i) = max(data.time{i});
end

% automatically determine the latency window which is possible in all trials
minperlength = [max(begsamplatency) min(endsamplatency)];
maxperlength = [min(begsamplatency) max(endsamplatency)];
maxtrllength = round((max(endsamplatency)-min(begsamplatency))*data.fsample) + 1;		% in samples
abstimvec    = ([1:maxtrllength] + min(data.offset) -1)./data.fsample;			      	% in seconds

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
        dat = preproc_baselinecorrect(dat);
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
        dat = preproc_baselinecorrect(dat);
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
cfg.version.id = '$Id: timelockanalysis.m,v 1.60 2009/03/23 21:21:16 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output 
timelock.cfg = cfg;

