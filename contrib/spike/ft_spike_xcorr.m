function [Xcorr] = ft_spike_xcorr(cfg,spike)

% FT_SPIKE_XCORR computes the cross-correlation histogram and shift predictor.
%
% Use as
%   [xcorr] = ft_spike_xcorr(cfg,data)
%
% The input SPIKE should be organised as the spike or the raw datatype, obtained from
% FT_SPIKE_MAKETRIALS or FT_PREPROCESSING (in that case, conversion is done
% within the function)
%
% Configurations options for xcorr general:
%   cfg.maxlag           = number in seconds, indicating the maximum lag for the
%                          cross-correlation function in sec (default = 0.1 sec).
%   cfg.biased           = 'yes' or 'no' (default). If 'no', we scale the
%                          cross-correlogram by M/(M-abs(lags)), where M = 2*N -1 with N
%                          the length of the data segment. Although this scaling reduces the
%                          bias, it can give a higher variance which might be more
%                          problematic in some cases.
%   cfg.shiftpredictor   = 'no' (default) or 'yes'. The shift-predictor is
%   calculated from
%                          channel x in every trial to the channel y in the previous
%                          trial. If two channels are independent, then the shift
%                          predictor should give the same correlogram as the raw
%                          correlogram calculated from the same trials.
%                          Typically, the shift predictor is subtracted from the
%                          correlogram.
%   cfg.outputunit       = 'proportion' (value in each bin indicates proportion of occurence)
%                          'center' (values are scaled to center value which is set to 1)
%                          'raw' (default) unnormalized crosscorrelogram.
%
%   cfg.binsize          = [binsize] in sec (default = 0.001 sec).
%
% Configuration options of channel selection & checking, latency and trial selection:
%   cfg.channelcmb       = Mx2 cell-array with selection of channel pairs (default = {'all' 'all'}),
%                          see FT_CHANNELCOMBINATION for details
%   cfg.latency          = [begin end] in seconds, 'max' (default), 'min', 'prestim'(t<=0), or
%                          'poststim' (t>=0).%
%   cfg.vartriallen      = 'yes' (default) or 'no'.
%                          If 'yes' - accept variable trial lengths and use all available trials
%                          and the samples in every trial.
%                          If 'no'  - only select those trials that fully cover the window as
%                          specified by cfg.latency and discard those trials that do not.
%   cfg.trials           = numeric selection of trials (default = 'all')
%   cfg.keeptrials       = 'yes' or 'no' (default)
%
% A peak at a negative lag for X.xcorr(:,chan1,chan2) means that chan1 is
% leading chan2. Thus, a negative lag represents a spike in the second
% dimension of X.xcorr before the channel in the third dimension of
% X.xcorr.
%
% Variable trial length is controlled by the option cfg.vartriallen. If
% cfg.vartriallen = 'yes', all trials are selected that have a minimum
% overlap with the latency window of cfg.maxlag. However, the shift
% predictor calculation demands that following trials have the same amount
% of data, otherwise, it does not control for rate non-stationarities. If
% cfg.vartriallen = 'yes', all trials should fall in the latency window,
% otherwise we do not compute the shift predictor.
%
% Output:
%   xcorr.xcorr            = (2*nlags+1)-by-nchans-by-nchans cross correlation histogram
%   xcorr.lags             = (2*nlags + 1) vector with lags in seconds.
%   xcorr.shiftpredictor   = (2*nlags+1)-by-nchans-by-nchans shift predictor.
%   xcorr.dimord
%   xcorr.trial            = (2*nlags + 1)-by-ntrials-by-nchans-by-nchans with single
%                            trials.
%   xcorr.label            = corresponding labels to channels in xcorr.xcorr
%   xcorr.cfg              = configurations used in this function.

% Copyright (C) 2010, Martin Vinck
%
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% check input spike structure
spike = ft_checkdata(spike,'datatype', 'spike', 'feedback', 'yes');

% get the default options
cfg.trials         = ft_getopt(cfg,'trials', 'all');
cfg.latency        = ft_getopt(cfg,'latency','maxperiod');
cfg.keeptrials     = ft_getopt(cfg,'keeptrials', 'yes');
cfg.shiftpredictor = ft_getopt(cfg,'shiftpredictor', 'no');
cfg.channelcmb     = ft_getopt(cfg,'channelcmb', 'all');
cfg.vartriallen    = ft_getopt(cfg,'vartriallen', 'no');
cfg.biased         = ft_getopt(cfg,'biased', 'no');
cfg.maxlag         = ft_getopt(cfg,'maxlag', 0.01);
cfg.binsize        = ft_getopt(cfg,'binsize', 0.001);
cfg.outputunit     = ft_getopt(cfg,'outputunit', 'proportion');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'latency', {'char', 'doublevector'});
cfg = ft_checkopt(cfg,'trials', {'char', 'doublevector', 'logical'}); 
cfg = ft_checkopt(cfg,'keeptrials', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'shiftpredictor', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'channelcmb', {'char', 'cell'});
cfg = ft_checkopt(cfg,'vartriallen', 'char', {'no', 'yes'});
cfg = ft_checkopt(cfg,'biased', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'maxlag', 'double');
cfg = ft_checkopt(cfg,'binsize', 'double');
cfg = ft_checkopt(cfg,'outputunit', 'char', {'proportion', 'center', 'raw'});

doShiftPredictor  = strcmp(cfg.shiftpredictor,'yes'); % shift predictor
doBiased          = strcmp(cfg.biased,'yes'); % debiasing

% determine the corresponding indices of the requested channel combinations
cfg.channelcmb = ft_channelcombination(cfg.channelcmb, spike.label,true);
cmbindx = zeros(size(cfg.channelcmb));
for k=1:size(cfg.channelcmb,1)
  cmbindx(k,1) = strmatch(cfg.channelcmb(k,1), spike.label, 'exact');
  cmbindx(k,2) = strmatch(cfg.channelcmb(k,2), spike.label, 'exact');
end
nCmbs 	   = size(cmbindx,1);
chansel    = unique(cmbindx(:)); % get the unique channels
nChans     = length(chansel);
if isempty(chansel),
  error('MATLAB:ft_spike_xcorr:cfg:channelcombination',...
    'No channel was selected')
end

% get the number of trials or change SPIKE according to cfg.trials
nTrials = size(spike.trialtime,1);
if  strcmp(cfg.trials,'all')
  cfg.trials = 1:nTrials;
elseif islogical(cfg.trials)
  cfg.trials = find(cfg.trials);
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>nTrials, error('MATLAB:ft_spike_xcorr:cfg:trials:maxExceeded',...
    'maximum trial number in cfg.trials should not exceed length of DATA.trial')
end
if isempty(cfg.trials),
  errors('MATLAB:fieldtrip:spike_xcorr:cfg:trials','No trials were selected in cfg.trials');
end
nTrials = length(cfg.trials);

% get the latencies
begTrialLatency = spike.trialtime(cfg.trials,1);
endTrialLatency = spike.trialtime(cfg.trials,2);
trialDur = endTrialLatency - begTrialLatency;

% select the latencies
if strcmp(cfg.latency,'minperiod')
  cfg.latency = [max(begTrialLatency) min(endTrialLatency)];
elseif strcmp(cfg.latency,'maxperiod')
  cfg.latency = [min(begTrialLatency) max(endTrialLatency)];
elseif strcmp(cfg.latency,'prestim')
  cfg.latency = [min(begTrialLatency) 0];
elseif strcmp(cfg.latency,'poststim')
  cfg.latency = [0 max(endTrialLatency)];
elseif ~isrealvec(cfg.latency)||length(cfg.latency)~=2
  error('MATLAB:ft_spike_xcorr:cfg:latency',...
    'cfg.latency should be "max", "min", "prestim", "poststim" or 1-by-2 numerical vector');
end
if cfg.latency(1)>cfg.latency(2), error('MATLAB:ft_spike_xcorr:incorrectLatencyWindow',...
    'cfg.latency should be a vector in ascending order, i.e., cfg.latency(2)>cfg.latency(1)');
end

% check whether the time window fits with the data
if (cfg.latency(1) < min(begTrialLatency)), cfg.latency(1) = min(begTrialLatency);
  warning('MATLAB:ft_spike_xcorr:incorrectLatencyWindow',...
    'Correcting begin latency of averaging window');
end
if (cfg.latency(2) > max(endTrialLatency)), cfg.latency(2) = max(endTrialLatency);
  warning('MATLAB:ft_spike_xcorr:incorrectLatencyWindow',...
    'Correcting begin latency of averaging window');
end

% get the maximum number of lags in samples
fsample = (1/cfg.binsize);
if cfg.maxlag>(cfg.latency(2)-cfg.latency(1))
  warning('MATLAB:ft_spike_xcorr:incorrectLatencyWindow',...
    'Correcting cfg.maxlag since it exceeds latency period');
  cfg.maxlag = (cfg.latency(2) - cfg.latency(1));
end
nLags        = round(cfg.maxlag*fsample);
lags         = (-nLags:nLags)/fsample; % create the lag axis

% check which trials will be used based on the latency
fullDur       = trialDur>=(cfg.maxlag); % only trials which are larger than maximum lag
overlaps      = endTrialLatency>(cfg.latency(1)+cfg.maxlag) & begTrialLatency<(cfg.latency(2)-cfg.maxlag);
hasWindow     = ones(nTrials,1);
if strcmp(cfg.vartriallen,'no') % only select trials that fully cover our latency window
  startsLater    = single(begTrialLatency>single(cfg.latency(1)));
  endsEarlier    = single(endTrialLatency<single(cfg.latency(2)));
  hasWindow      = ~(startsLater | endsEarlier); % check this in all other funcs
end
trialSel           = fullDur(:) & overlaps(:) & hasWindow(:);
cfg.trials         = cfg.trials(trialSel);
if isempty(cfg.trials),
  warning('MATLAB:ft_spike_xcorr:cfg:trials','No trials were selected');
end
if length(cfg.trials)<2&&doShiftPredictor
  warning('MATLAB:ft_spike_xcorr:cfg:shiftpredictor:trials',...
    'Shift predictor can only be calculated with more than 1 selected trial')
  doShiftPredictor = 0;
end

% if we allow variable trial lengths
if strcmp(cfg.vartriallen,'yes') && doShiftPredictor && nTrials~=length(cfg.trials)
  warning('MATLAB:ft_spike_xcorr:cfg:vartriallen:shiftpredictor:latency',...
    '%s\n%s\n%s', 'Variable trial length and shift predictor are only possible when all trials',...
    'contain the full window period, ignoring request to compute shift predictor',...
    'please adjust the latency so it includes all trials or set vartriallen to no');
  doShiftPredictor = 0;
end
nTrials   = length(cfg.trials); % only now reset nTrials

% preallocate the sum and the single trials and the shift predictor
s            = zeros(2*nLags + 1,nChans,nChans);
keepTrials   = strcmp(cfg.keeptrials,'yes');
if keepTrials,       singleTrials = zeros(2*nLags+1,nTrials,nChans,nChans);      end
if doShiftPredictor, shiftSum     = zeros(2*nLags + 1,nChans,nChans);            end
for iTrial = 1:nTrials
  origTrial = cfg.trials(iTrial);

  for iCmb = 1:nCmbs
    % take only the times that are in this trial
    indx  = cmbindx(iCmb,:);
    inTrial1 = spike.trial{indx(1)}==origTrial;
    inTrial2 = spike.trial{indx(2)}==origTrial;
    ts1 = sort(spike.time{indx(1)}(inTrial1));
    ts2 = sort(spike.time{indx(2)}(inTrial2));

    % compute the xcorr if both are non-empty
    if ~isempty(ts1) && ~isempty(ts2)
      if indx(1)<=indx(2)
        [x]   = spike_crossx(ts1(:),ts2(:),cfg.binsize,nLags*2+1);
      else
        [x]   = spike_crossx(ts2(:),ts1(:),cfg.binsize,nLags*2+1);
      end

      % sum the xcorr
      s(:,indx(1),indx(2)) = s(:,indx(1),indx(2)) + x(:);
      s(:,indx(2),indx(1)) = s(:,indx(2),indx(1)) + flipud(x(:));

      % store individual trials if requested
      if keepTrials
        singleTrials(:,iTrial,indx(1),indx(2)) = x(:);
        singleTrials(:,iTrial,indx(2),indx(1)) = flipud(x(:));
      end
    end

    % compute the shift predictor
    if doShiftPredictor && iTrial>1

      % symmetric, get x21 to x12 and x22 to x11
      inTrial1_old = spike.trial{indx(1)}==cfg.trials(iTrial-1);
      ts1_old      = sort(spike.time{indx(1)}(inTrial1_old));
      inTrial2_old = spike.trial{indx(2)}==cfg.trials(iTrial-1);
      ts2_old      = sort(spike.time{indx(2)}(inTrial2_old));

      % compute both combinations
      for k = 1:2

        % take chan from this and previous channel
        if k==1
          A = ts1;
          B = ts2_old;
        else
          A = ts1_old;
          B = ts2;
        end
        if ~isempty(A) && ~isempty(B),

          ts1(:)
          ts2(:)
          if indx(1)<=indx(2)            
            [x]   = spike_crossx(A(:),B(:),cfg.binsize,nLags*2+1);
          else
            [x]   = spike_crossx(A(:),B(:),cfg.binsize,nLags*2+1);
          end
          % compute the sum
          shiftSum(:,indx(1),indx(2)) =  shiftSum(:,indx(1),indx(2)) + x(:);
          shiftSum(:,indx(2),indx(1)) =  shiftSum(:,indx(2),indx(1)) + flipud(x(:));
        end
      end
    end % symmetric shift predictor loop
  end % combinations
end % trial loop

% multiply the shift sum by a factor so it has the same scale: note it is not raw anymore
dofShiftPred = 2*(nTrials-1);
if doShiftPredictor, shiftSum = shiftSum*nTrials/dofShiftPred; end

switch cfg.outputunit
  case 'proportion'
    if doShiftPredictor
      sm       = sum(shiftSum);
      shiftSum = shiftSum./sm(ones(nLags*2+1,1),:,:);
    end
    
    sm       = sum(s);
    s        = s./sm(ones(nLags*2+1,1),:,:);
  case 'center'
    if doShiftPredictor
      center  = shiftSum(nLags+1,:,:);
      shiftSum = shiftSum./center(ones(nLags*2+1,1),:,:);
    end
    center   = sum(s);
    s        = s./center(ones(nLags*2+1,1),:,:);
end

% return the results
Xcorr.xcorr          = s;
Xcorr.lags           = lags;
if doShiftPredictor
  Xcorr.shiftpredictor = shiftSum;
end
if keepTrials
  Xcorr.trial = singleTrials;
end
Xcorr.label       = spike.label(chansel);
Xcorr.labelcmb    = cfg.channelcmb;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous spike
ft_postamble history Xcorr

