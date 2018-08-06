function [sdf, sdfdata] = ft_spikedensity(cfg, data)

% FT_SPIKEDENSITY computes the spike density function of the spike trains by
% convolving the data with a window.
%
% Use as
%   [sdf]          = ft_spike_density(cfg, data)
%   [sdf, sdfdata] = ft_spike_density(cfg, data)
% 
% If you specify one output argument, only the average and variance of spike density
% function across trials will be computed and individual trials are not kept. See
% cfg.winfunc below for more information on the smoothing kernel to use.
%
% Inputs:
%   DATA should be organised in a RAW structure with binary spike
%   representations obtained from FT_APPENDSPIKE or FT_CHECKDATA, or
%   a SPIKE structure.
%
% Configurations:
%   cfg.timwin         = [begin end], time of the smoothing kernel (default = [-0.05 0.05])
%                        If cfg.winfunc = @alphawin, cfg.timwin(1) will be
%                        set to 0. Hence, it is possible to use asymmetric
%                        kernels. 
%   cfg.outputunit     = 'rate' (default) or 'spikecount'. This determines the physical unit
%                        of our spikedensityfunction, either in firing rate or in spikecount.
%   cfg.winfunc        = (a) string or function handle, type of window to convolve with (def = 'gauss').
%                        - 'gauss' (default)
%                        - 'alphawin', given by win = x*exp(-x/timeconstant)
%                        - For standard window functions in the signal processing toolbox see
%                          WINDOW.
%                        (b) vector of length nSamples, used directly as window
%   cfg.winfuncopt     = options that go with cfg.winfunc
%                        For cfg.winfunc = 'alpha': the timeconstant in seconds (default = 0.005s)
%                        For cfg.winfunc = 'gauss': the standard deviation in seconds (default =
%                                         1/4 of window duration in seconds)
%                        For cfg.winfunc = 'wname' with 'wname' any standard window function
%                                          see window opts in that function and add as cell array
%                        If cfg.winfunctopt = [], default opts are taken.
%   cfg.latency        = [begin end] in seconds, 'maxperiod' (default), 'minperiod',
%                        'prestim'(t>=0), or 'poststim' (t>=0).
%   cfg.vartriallen    = 'yes' (default) or 'no'.
%                        'yes' - accept variable trial lengths and use all available trials
%                         and the samples in every trial. Missing values will be ignored in
%                         the computation of the average and the variance.
%                        'no'  - only select those trials that fully cover the window as
%                         specified by cfg.latency.
%   cfg.spikechannel   = cell-array ,see FT_CHANNELSELECTION for details
%   cfg.trials         =  numeric or logical selection of trials (default = 'all')
%   cfg.keeptrials     = 'yes' or 'no' (default). If 'yes', we store the trials in a matrix
%                         in the output SDF as well
%   cfg.fsample        = additional user input that can be used when input
%                        is a SPIKE structure, in that case a continuous
%                        representation is created using cfg.fsample
%                        (default = 1000)
%
% The SDF output is a data structure similar to the TIMELOCK structure from FT_TIMELOCKANALYSIS.
% For subsequent processing you can use for example
%   FT_TIMELOCKSTATISTICS                Compute statistics on SDF
%   FT_SPIKE_PLOT_RASTER                 Plot together with the raster plots
%   FT_SINGLEPLOTER and FT_MULTIPLOTER   Plot spike-density alone
%
% The SDFDATA output is a data structure similar to DATA type structure from FT_PREPROCESSING.
% For subsequent processing you can use for example
%   FT_TIMELOCKANALYSIS                  Compute time-locked average and variance
%   FT_FREQANALYSIS                      Compute frequency and time-ferquency spectrum.

% Copyright (C) 2010-2012, Martin Vinck
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
ft_preamble provenance data
ft_preamble trackconfig

% get the default options
if isfield(cfg,'trials') && isempty(cfg.trials), error('no trials were selected'); end % empty should result in error, not in default
cfg.outputunit   = ft_getopt(cfg,'outputunit','rate');
cfg.timwin       = ft_getopt(cfg,'timwin',[-0.05 0.05]);
cfg.trials       = ft_getopt(cfg,'trials', 'all');
cfg.latency      = ft_getopt(cfg,'latency','maxperiod');
cfg.spikechannel = ft_getopt(cfg,'spikechannel', 'all');
cfg.vartriallen  = ft_getopt(cfg,'vartriallen', 'yes');
cfg.keeptrials   = ft_getopt(cfg,'keeptrials', 'yes');
cfg.winfunc      = ft_getopt(cfg,'winfunc', 'gauss');
cfg.winfuncopt   = ft_getopt(cfg,'winfuncopt', []);
cfg.fsample      = ft_getopt(cfg,'fsample', 1000);

% ensure that the options are valid
cfg = ft_checkopt(cfg, 'outputunit','char', {'rate', 'spikecount'});
cfg = ft_checkopt(cfg, 'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg, 'latency', {'char', 'ascendingdoublebivector'});
cfg = ft_checkopt(cfg, 'trials', {'char', 'doublevector', 'logical'}); 
cfg = ft_checkopt(cfg, 'vartriallen', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg, 'keeptrials', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg, 'timwin', 'ascendingdoublebivector');
cfg = ft_checkopt(cfg, 'winfunc', {'char', 'function_handle', 'doublevector'});
cfg = ft_checkopt(cfg, 'winfuncopt', {'cell', 'double', 'empty'});
cfg = ft_checkopt(cfg, 'fsample', 'double');
if strcmp(class(cfg.winfunc), 'function_handle'), cfg.winfunc = func2str(cfg.winfunc); end

cfg = ft_checkconfig(cfg, 'allowed', {'outputunit', 'spikechannel', 'latency', 'trials', 'vartriallen', 'keeptrials', 'timwin', 'winfunc', 'winfuncopt', 'fsample'});

% check input data structure
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes', 'fsample', cfg.fsample);

[spikechannel] = detectspikechan(data);
if strcmp(cfg.spikechannel, 'all'), 
  cfg.spikechannel = spikechannel; 
else
  cfg.spikechannel = ft_channelselection(cfg.spikechannel, data.label);  
  if ~all(ismember(cfg.spikechannel,spikechannel)), warning('some selected spike channels no not appear spike channels'); end        
end
spikesel    = match_str(data.label, cfg.spikechannel);
nUnits      = length(spikesel); % number of spike channels
if nUnits==0, error('no spikechannel selected by means of cfg.spikechannel'); end

% get the number of trials or change DATA according to cfg.trials
if  strcmp(cfg.trials,'all')
  cfg.trials = 1:length(data.trial);
elseif islogical(cfg.trials)
  cfg.trials = find(cfg.trials);
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>length(data.trial),error('maximum trial number in cfg.trials should not exceed length of data.trial'), end
if isempty(cfg.trials), error('no trials were selected in cfg.trials'); end

% determine the duration of each trial
begTrialLatency = cellfun(@min,data.time(cfg.trials));
endTrialLatency = cellfun(@max,data.time(cfg.trials));

% select the latencies
if strcmp(cfg.latency,'minperiod')
  cfg.latency = [max(begTrialLatency) min(endTrialLatency)];
elseif strcmp(cfg.latency,'maxperiod')
  cfg.latency = [min(begTrialLatency) max(endTrialLatency)];
elseif strcmp(cfg.latency,'prestim')
  cfg.latency = [min(begTrialLatency) 0];
elseif strcmp(cfg.latency,'poststim')
  cfg.latency = [0 max(endTrialLatency)];
end
if (cfg.latency(1) < min(begTrialLatency)), cfg.latency(1) = min(begTrialLatency);
  warning('correcting begin latency of averaging window');
end
if (cfg.latency(2) > max(endTrialLatency)), cfg.latency(2) = max(endTrialLatency);
  warning('correcting end latency of averaging window');
end

% start processing the window information
if strcmp(cfg.winfunc,'alphawin') % now force start of window to be positive.
  warning('cfg.timwin(1) should be a non-negative number if cfg.winfunc = @alphawin, correcting')
  cfg.timwin(1) = 0;
end
if cfg.timwin(1)>0 || cfg.timwin(2)<0, error('Please specify cfg.timwin(1)<=0 and cfg.timwin(2)>=0'); end

% construct the window and the time axis of the window
fsample       = data.fsample; 
sampleTime    = 1/fsample;
nLeftSamples  = round(-cfg.timwin(1)/sampleTime);
nRightSamples = round(cfg.timwin(2)/sampleTime);
winTime       = -nLeftSamples*sampleTime : sampleTime : nRightSamples*sampleTime; % this is uneven if cfg.timwin symmetric
cfg.timwin    = [winTime(1) winTime(end)];
nSamplesWin   = length(winTime);
if nSamplesWin==1, warning('Number of samples in selected window is exactly one, so no smoothing applied'); end

% construct the window
if ~iscell(cfg.winfuncopt), cfg.winfuncopt = {cfg.winfuncopt}; end
if strcmp(cfg.winfunc,'gauss')
  if  isempty(cfg.winfuncopt{1}), cfg.winfuncopt{1} = 0.25*diff(cfg.timwin); end
  win = exp(-(winTime.^2)/(2*cfg.winfuncopt{1}^2)); % here we could compute the optimal SD
elseif strcmp(cfg.winfunc,'alphawin')
  if isempty(cfg.winfuncopt{1}),  cfg.winfuncopt{1} = 0.005; end
  win = winTime.*exp(-winTime/cfg.winfuncopt{1});
elseif ischar(cfg.winfunc) 
  if isempty(cfg.winfuncopt{1})
    win = feval(cfg.winfunc,nSamplesWin);
  else
    win = feval(cfg.winfunc,nSamplesWin, cfg.winfuncopt{:});
  end
else % must be a double vector then
  if length(cfg.winfunc)~=nSamplesWin, error('Length of cfg.winfunc vector should be %d', nSamplesWin); end
end
win    = win(:).'./sum(win);    % normalize the window to 1
winDur = max(winTime) - min(winTime); % duration of the window

% create the time axis for the spike density
time           = cfg.latency(1):(1/fsample):cfg.latency(2);
cfg.latency(2) = time(end);   % this is the used latency that should be stored in cfg again
maxNumSamples  = length(time);

% check which trials will be used based on the latency
overlaps      = endTrialLatency>=(cfg.latency(1)+winDur) & begTrialLatency<=(cfg.latency(2)-winDur);
if strcmp(cfg.vartriallen,'no') % only select trials that fully cover our latency window
  hasWindow       = false(length(begTrialLatency),1);
  for iTrial      = 1:length(begTrialLatency)
    timeTrial     = data.time{cfg.trials(iTrial)};
    nSamplesTrial = length(timeTrial);
    hasLaterStart = (begTrialLatency(iTrial)-cfg.latency(1))>(0.5*sampleTime);    
    hasEarlierEnd = (cfg.latency(2)-endTrialLatency(iTrial))>(0.5*sampleTime);
    hasWindow(iTrial) = ~hasLaterStart & ~hasEarlierEnd & nSamplesTrial>=maxNumSamples;
  end
else
  hasWindow      = true(length(cfg.trials),1);  % in case vartriallen = "yes"
end
trialSel          = overlaps(:) & hasWindow(:); % trials from cfg.trials we select further
cfg.trials        = cfg.trials(trialSel);       % cut down further on cfg.trials
begTrialLatency   = begTrialLatency(trialSel);  % on this variable as well
endTrialLatency   = endTrialLatency(trialSel);  % on this variable as well
nTrials           = length(cfg.trials);         % the actual number of trials we will use
if isempty(cfg.trials),warning('no trials were selected, please check cfg.trials'), end

% calculates the samples we are shifted wrt latency(1)
samplesShift      = zeros(1,nTrials);
sel               = (begTrialLatency-cfg.latency(1))>(0.5*sampleTime); % otherwise 0 samples to be padded
samplesShift(sel) = round(fsample*(begTrialLatency(sel)-cfg.latency(1)));

% preallocate the sum, squared sum and degrees of freedom
[s,ss]   = deal(NaN(nUnits, maxNumSamples)); % sum and sum of squares
dof      = zeros(nUnits, length(s));

% preallocate, depending on whether nargout is 1 or 2
if (strcmp(cfg.keeptrials,'yes')), singleTrials = zeros(nTrials,nUnits,size(s,2)); end
if nargout==2, [sdfdata.trial(1:nTrials) sdfdata.time(1:nTrials)] = deal({[]}); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           compute the spike density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iTrial = 1:nTrials
  origTrial  = cfg.trials(iTrial);   % this is the original trial we use for DATA input
  timeAxis   = data.time{origTrial}; % get the time axis for this trial
  begSample  = nearest(timeAxis, cfg.latency(1));
  sampleSel  = begSample : (begSample + maxNumSamples - samplesShift(iTrial) - 1); 
  sampleSel(sampleSel>length(timeAxis)) = []; % delete the indices that should not be there
  nSamples   = length(sampleSel);
  trialTime  = timeAxis(sampleSel);  % select the relevant portion of time
  
  % handle every unit separately
  for iUnit = 1:nUnits
    unitIndx = spikesel(iUnit); % index in data.label
    dat      = data.trial{origTrial}(unitIndx,sampleSel); % get the data
    if any(dat)
      y = conv(full(dat),win);  % use convolution to get the raw spike density
    else
      y = zeros(1,nSamples + nSamplesWin - 1); % no spikes; no convolution needed
    end
    
    if strcmp(cfg.outputunit, 'rate')
      y = y*fsample;  % normalize to the sampling rate, to get it in Hz.
    else
      y = y*nSamplesWin;   % now maximum becomes N (length window)
    end
    
    % restrict to the relevant portion of output conv
    y = y(nLeftSamples+1 : end-nRightSamples); % delete additional points we get with conv
    y([1:nLeftSamples end-nRightSamples+1:end]) = NaN; % assign NaN at borders
    
    % write a raw data structure
    if nargout==2
      sl = ~isnan(y);
      sdfdata.trial{iTrial}(iUnit,:) = y(sl);
      sdfdata.time{iTrial}(1,:)      = trialTime(sl); % write back original time axis
    end
    
    % pad with nans if there's variable trial length
    dofsel = ~isnan(y);%true(1,length(y));
    if strcmp(cfg.vartriallen,'yes')
      padLeft  = zeros(1, samplesShift(iTrial));
      padRight = zeros(1,(maxNumSamples - nSamples - samplesShift(iTrial)));
      ySingleTrial = [NaN(size(padLeft)) y NaN(size(padRight))];
      y        = [padLeft y padRight];
      dofsel   = logical([padLeft dofsel padRight]);
    else
      ySingleTrial = y;
    end
    
    % compute the sum and the sum of squares
    s(iUnit,:)        = nansum([s(iUnit,:);y]);      % compute the sum
    ss(iUnit,:)       = nansum([ss(iUnit,:);y.^2]);  % compute the squared sum
    
    % count the number of samples that went into the sum
    dof(iUnit,dofsel) = dof(iUnit,dofsel) + 1;
    
    % keep the single trials if requested
    if strcmp(cfg.keeptrials,'yes'), singleTrials(iTrial,iUnit,:) = ySingleTrial; end
  end
  
  % remove the trial from data in order to avoid buildup in memory
  data.trial{origTrial} = [];
  data.time{origTrial}  = [];
end

% give back a similar structure as timelockanalysis
sdf.avg                  = s ./ dof;
sdf.var                  = (ss - s.^2./dof)./(dof-1); 
sdf.dof                  = dof;
sdf.time                 = time;
sdf.label(1:nUnits)      = data.label(spikesel);
if (strcmp(cfg.keeptrials,'yes'))
  sdf.trial = singleTrials;
  sdf.dimord = 'rpt_chan_time';
else
  sdf.dimord = 'chan_time';
end
  
% create a new structure that is a standard raw data spike structure itself, this is returned as second output argument
if nargout==2
  sdfdata.fsample                   = fsample;
  sdfdata.label(1:nUnits)           = data.label(spikesel);
  try, sdfdata.hdr                  = data.hdr; end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble previous data
ft_postamble provenance sfd
ft_postamble history    sdf
if nargout==2
  ft_postamble history sdfdata
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spikelabel, eeglabel] = detectspikechan(data)

% autodetect the spike channels
ntrial = length(data.trial);
nchans  = length(data.label);
spikechan = zeros(nchans,1);
for i=1:ntrial
  for j=1:nchans
    hasAllInts    = all(isnan(data.trial{i}(j,:)) | data.trial{i}(j,:) == round(data.trial{i}(j,:)));
    hasAllPosInts = all(isnan(data.trial{i}(j,:)) | data.trial{i}(j,:)>=0);
    spikechan(j)  = spikechan(j) + double(hasAllInts & hasAllPosInts);
  end
end
spikechan = (spikechan==ntrial);

spikelabel = data.label(spikechan);
eeglabel   = data.label(~spikechan);

