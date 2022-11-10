function [stat] = ft_spike_xcorr(cfg, spike)

% FT_SPIKE_XCORR computes the cross-correlation histogram and shift predictor.
%
% Use as
%   [stat] = ft_spike_xcorr(cfg, data)
%
% The input SPIKE should be organised as the spike or the raw datatype, obtained from
% FT_SPIKE_MAKETRIALS or FT_PREPROCESSING (in that case, conversion is done
% within the function). A mex file is located in /contrib/spike/private
% which will be automatically mexed.
%
% Configurations options for xcorr general:
%   cfg.maxlag           = number in seconds, indicating the maximum lag for the
%                          cross-correlation function in sec (default = 0.1 sec).
%   cfg.debias           = 'yes' (default) or 'no'. If 'yes', we scale the
%                          cross-correlogram by M/(M-abs(lags)), where M = 2*N -1 with N
%                          the length of the data segment.
%   cfg.method           = 'xcorr' or 'shiftpredictor'. If 'shiftpredictor'
%                           we do not compute the normal cross-correlation
%                           but shuffle the subsequent trials.
%                           If two channels are independent, then the shift
%                           predictor should give the same correlogram as the raw
%                           correlogram calculated from the same trials.
%                           Typically, the shift predictor is subtracted from the
%                           correlogram.
%   cfg.outputunit       = - 'proportion' (value in each bin indicates proportion of occurence)
%                          - 'center' (values are scaled to center value which is set to 1)
%                          - 'raw' (default) unnormalized crosscorrelogram.
%   cfg.binsize          = [binsize] in sec (default = 0.001 sec).
%   cfg.channelcmb       = Mx2 cell-array with selection of channel pairs (default = {'all' 'all'}),
%                          see FT_CHANNELCOMBINATION for details
%   cfg.latency          = [begin end] in seconds, 'max' (default), 'min', 'prestim'(t<=0), or
%                          'poststim' (t>=0).%
%   cfg.vartriallen      = 'yes' (default) or 'no'.
%                          If 'yes' - accept variable trial lengths and use all available trials
%                          and the samples in every trial.
%                          If 'no'  - only select those trials that fully cover the window as
%                          specified by cfg.latency and discard those trials that do not.
%                          if cfg.method = 'yes', then cfg.vartriallen
%                          should be 'no' (otherwise, fewer coincidences
%                          will occur because of non-overlapping windows)
%   cfg.trials           = numeric selection of trials (default = 'all')
%   cfg.keeptrials       = 'yes' or 'no' (default)
%
% A peak at a negative lag for stat.xcorr(chan1,chan2,:) means that chan1 is leading
% chan2. Thus, a negative lag represents a spike in the second dimension of
% stat.xcorr before the channel in the third dimension of stat.stat.
%
% Variable trial length is controlled by the option cfg.vartriallen. If it is
% specified as cfg.vartriallen='yes', all trials are selected that have a minimum
% overlap with the latency window of cfg.maxlag. However, the shift predictor
% calculation demands that following trials have the same amount of data, otherwise,
% it does not control for rate non-stationarities. If cfg.vartriallen = 'yes', all
% trials should fall in the latency window, otherwise we do not compute the shift
% predictor.
%
% Output:
%    stat.xcorr            = nchans-by-nchans-by-(2*nlags+1) cross correlation histogram with dimord 'chan_chan_time'
%   or
%    stat.shiftpredictor   = nchans-by-nchans-by-(2*nlags+1) shift predictor with dimord 'chan_chan_time'
% and
%    stat.lags             = (2*nlags + 1) vector with lags in seconds.
%    stat.trial            = ntrials-by-nchans-by-nchans-by-(2*nlags + 1) with single trials and dimord 'rpt_chan_chan_time'
%    stat.label            = corresponding labels to channels in stat.xcorr
%    stat.cfg              = configurations used in this function

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
ft_preamble provenance spike


% check input spike structure
spike = ft_checkdata(spike, 'datatype', 'spike', 'feedback', 'yes');

% get the default options
cfg.trials         = ft_getopt(cfg, 'trials', 'all');
cfg.latency        = ft_getopt(cfg, 'latency','maxperiod');
cfg.keeptrials     = ft_getopt(cfg, 'keeptrials', 'no');
cfg.method         = ft_getopt(cfg, 'method', 'xcorr');
cfg.channelcmb     = ft_getopt(cfg, 'channelcmb', 'all');
cfg.vartriallen    = ft_getopt(cfg, 'vartriallen', 'yes');
cfg.debias         = ft_getopt(cfg, 'debias', 'yes');
cfg.maxlag         = ft_getopt(cfg, 'maxlag', 0.01);
cfg.binsize        = ft_getopt(cfg, 'binsize', 0.001);
cfg.outputunit     = ft_getopt(cfg, 'outputunit', 'proportion');

% ensure that the options are valid
cfg = ft_checkopt(cfg, 'latency', {'char', 'ascendingdoublebivector'});
cfg = ft_checkopt(cfg, 'trials', {'char', 'doublevector', 'logical'});
cfg = ft_checkopt(cfg, 'keeptrials', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg, 'method', 'char', {'xcorr', 'shiftpredictor'});
cfg = ft_checkopt(cfg, 'channelcmb', {'char', 'cell'});
cfg = ft_checkopt(cfg, 'vartriallen', 'char', {'no', 'yes'});
cfg = ft_checkopt(cfg, 'debias', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg, 'maxlag', 'doublescalar');
cfg = ft_checkopt(cfg, 'binsize', 'doublescalar');
cfg = ft_checkopt(cfg, 'outputunit', 'char', {'proportion', 'center', 'raw'});

cfg = ft_checkconfig(cfg, 'allowed', {'latency', 'trials', 'keeptrials', 'method', 'channelcmb', 'vartriallen', 'debias', 'maxlag', 'binsize', 'outputunit'});

doShiftPredictor  = strcmp(cfg.method, 'shiftpredictor'); % shift predictor

% determine the corresponding indices of the requested channel combinations
cfg.channelcmb = ft_channelcombination(cfg.channelcmb, spike.label(:), true);
cmbindx        = zeros(size(cfg.channelcmb));
for k=1:size(cfg.channelcmb,1)
  cmbindx(k,1) = strmatch(cfg.channelcmb(k,1), spike.label, 'exact');
  cmbindx(k,2) = strmatch(cfg.channelcmb(k,2), spike.label, 'exact');
end
nCmbs 	   = size(cmbindx,1);
chansel    = unique(cmbindx(:)); % get the unique channels
nChans     = length(chansel);
if isempty(chansel), error('No channel was selected'); end

% get the number of trials or change SPIKE according to cfg.trials
nTrials = size(spike.trialtime,1);
if  strcmp(cfg.trials,'all')
  cfg.trials = 1:nTrials;
elseif islogical(cfg.trials) || all(cfg.trials==0 | cfg.trials==1)
  cfg.trials = find(cfg.trials);
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>nTrials, error('maximum trial number in cfg.trials should not exceed length of DATA.trial'); end
if isempty(cfg.trials), errors('No trials were selected in cfg.trials'); end
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
end

% check whether the time window fits with the data
if (cfg.latency(1) < min(begTrialLatency)), cfg.latency(1) = min(begTrialLatency);
  warning('Correcting begin latency of averaging window');
end
if (cfg.latency(2) > max(endTrialLatency)), cfg.latency(2) = max(endTrialLatency);
  warning('Correcting begin latency of averaging window');
end

% get the maximum number of lags in samples
fsample = (1/cfg.binsize);
if cfg.maxlag>(cfg.latency(2)-cfg.latency(1))
  warning('Correcting cfg.maxlag since it exceeds latency period');
  cfg.maxlag = (cfg.latency(2) - cfg.latency(1));
end
nLags        = round(cfg.maxlag*fsample);
lags         = (-nLags:nLags)/fsample; % create the lag axis

% check which trials will be used based on the latency
fullDur       = trialDur>=(cfg.maxlag); % only trials which are larger than maximum lag
overlaps      = endTrialLatency>(cfg.latency(1)+cfg.maxlag) & begTrialLatency<(cfg.latency(2)-cfg.maxlag);
hasWindow     = ones(nTrials,1);
if strcmp(cfg.vartriallen,'no') % only select trials that fully cover our latency window
  startsLater    = single(begTrialLatency) > (single(cfg.latency(1)) + 0.5*cfg.binsize);
  endsEarlier    = single(endTrialLatency) < (single(cfg.latency(2)) - 0.5*cfg.binsize);
  hasWindow      = ~(startsLater | endsEarlier); % check this in all other funcs
end
trialSel           = fullDur(:) & overlaps(:) & hasWindow(:);
cfg.trials         = cfg.trials(trialSel);
begTrialLatency    = begTrialLatency(trialSel);
endTrialLatency    = endTrialLatency(trialSel);

if isempty(cfg.trials), warning('No trials were selected'); end
if length(cfg.trials)<2&&doShiftPredictor
  error('Shift predictor can only be calculated with more than 1 selected trial.');
end

% if we allow variable trial lengths
if strcmp(cfg.vartriallen,'yes') && doShiftPredictor && ~strcmp(cfg.outputunit,'proportion')
  warning('using cfg.vartriallen = "yes" and shift predictor method: please use cfg.outputunit = "proportion"');
end
nTrials   = length(cfg.trials); % only now reset nTrials

% preallocate the sum and the single trials and the shift predictor
keepTrials   = strcmp(cfg.keeptrials,'yes');
s     = zeros(nChans,nChans,2*nLags);
if keepTrials
  warning('storing single trials for cross correlation is memory expensive, please check');
  singleTrials = zeros(nTrials,nChans,nChans,2*nLags);
end

if strcmp(cfg.method,'shiftpredictor'), singleTrials(1,:,:,:) = NaN; end

if ((cfg.latency(2)-cfg.latency(1))/cfg.maxlag)<2.5
  warning('selected latency will cause highly variable xcorr at borders');
end

ft_progress('init', 'text',     'Please wait...');
for iTrial = 1:nTrials
  origTrial = cfg.trials(iTrial);
  ft_progress(iTrial/nTrials, 'Processing trial %d from %d', iTrial, nTrials);
  
  for iCmb = 1:nCmbs
    % take only the times that are in this trial and latency window
    indx  = cmbindx(iCmb,:);
    inTrial1 = spike.trial{indx(1)}==origTrial;
    inTrial2 = spike.trial{indx(2)}==origTrial;
    inWindow1 = spike.time{indx(1)}>=cfg.latency(1) &  spike.time{indx(1)}<=cfg.latency(2);
    inWindow2 = spike.time{indx(2)}>=cfg.latency(1) &  spike.time{indx(2)}<=cfg.latency(2);
    ts1 = sort(spike.time{indx(1)}(inTrial1(:) & inWindow1(:)));
    ts2 = sort(spike.time{indx(2)}(inTrial2(:) & inWindow2(:)));
    
    switch cfg.method
      case 'xcorr'
        % compute the xcorr if both are non-empty
        if ~isempty(ts1) && ~isempty(ts2)
          if indx(1)<=indx(2)
            [x]   = spike_crossx_matlab(ts1(:),ts2(:),cfg.binsize,nLags*2+1);  % removed the mex file for now, until issue on windows is resolved
          else
            [x]   = spike_crossx_matlab(ts2(:),ts1(:),cfg.binsize,nLags*2+1); % removed the mex file for now, until issue on windows is resolved
          end
          
          % remove the center peaks from the auto-correlogram
          if indx(1)==indx(2), x(nLags:nLags+1) = 0; end
          
          if strcmp(cfg.debias,'yes')
            lags = (-nLags:nLags)*cfg.binsize;
            lags = (lags(2:end)+lags(1:end-1))/2;
            T    = nanmin([endTrialLatency(iTrial);cfg.latency(2)])-nanmax([begTrialLatency(iTrial);cfg.latency(1)]);
            sc = (T./(T-abs(lags(:))));
            sc = length(sc)*sc./sum(sc);
            x  = x(:).*sc(:);
          end
          % sum the xcorr
          s(indx(1),indx(2),:) = s(indx(1),indx(2),:) + shiftdim(x(:),-2);
          s(indx(2),indx(1),:) = s(indx(2),indx(1),:) + shiftdim(flipud(x(:)),-2);
          
          % store individual trials if requested
          if keepTrials
            singleTrials(iTrial,indx(1),indx(2),:) = shiftdim(x(:),-3);
            singleTrials(iTrial,indx(2),indx(1),:) = shiftdim(flipud(x(:)),-3);
          end
        end
      case 'shiftpredictor'
        if iTrial>1
          % symmetric, get x21 to x12 and x22 to x11
          inTrial1_old = spike.trial{indx(1)}==cfg.trials(iTrial-1);
          ts1_old      = sort(spike.time{indx(1)}(inTrial1_old(:) & inWindow1(:)));
          inTrial2_old = spike.trial{indx(2)}==cfg.trials(iTrial-1);
          ts2_old      = sort(spike.time{indx(2)}(inTrial2_old(:) & inWindow2(:)));
          
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
            if ~isempty(A) && ~isempty(B)
              if indx(1)<=indx(2)
                [x]   = spike_crossx_matlab(A(:),B(:),cfg.binsize,nLags*2+1);
              else
                [x]   = spike_crossx_matlab(B(:),A(:),cfg.binsize,nLags*2+1);
              end
              % remove the center peaks from the auto-correlogram
              if indx(1)==indx(2), x(nLags:nLags+1) = 0; end
              
              if strcmp(cfg.debias,'yes')
                lags = (-nLags:nLags)*cfg.binsize;
                lags = (lags(2:end)+lags(1:end-1))/2;
                T    = nanmin([endTrialLatency(iTrial);cfg.latency(2)])-nanmax([begTrialLatency(iTrial);cfg.latency(1)]);
                sc = (T./(T-abs(lags(:))));
                sc = length(sc)*sc./sum(sc);
                x  = x(:).*sc(:);
              end
              % compute the sum
              s(indx(1),indx(2),:) =  s(indx(1),indx(2),:) + shiftdim(x(:),-2);
              s(indx(2),indx(1),:) =  s(indx(2),indx(1),:) + shiftdim(flipud(x(:)),-2);
              if keepTrials
                singleTrials(iTrial,indx(1),indx(2),:) = shiftdim(x(:)/2,-3);
                singleTrials(iTrial,indx(2),indx(1),:) = shiftdim(flipud(x(:))/2,-3);
              end
            end
          end
        end
    end % symmetric shift predictor loop
  end % combinations
end % trial loop
ft_progress('close')
% multiply the shift sum by a factor so it has the same scale as raw
dofShiftPred = 2*(nTrials-1);
if doShiftPredictor, s = s*nTrials/dofShiftPred; end
switch cfg.outputunit
  case 'proportion'
    sm = repmat(nansum(s,3),[1 1 nLags*2]);
    s = s./sm;
    if keepTrials
      singleTrials       = singleTrials./repmat(shiftdim(sm,-1),[nTrials 1 1 1]);
    end
  case 'center'
    center   = repmat(s(:,:,nLags+1),[1 1 nLags*2]);
    s        = s./center;
    if keepTrials
      singleTrials       = singleTrials./repmat(shiftdim(center,-1),[nTrials 1 1 1]);
    end
end

% return the results
stat.(cfg.method)      = s;
lags = (-nLags:nLags)*cfg.binsize;
lags = (lags(2:end)+lags(1:end-1))/2;
stat.time              = lags;
stat.dimord            = 'chan_chan_time';
if keepTrials
  stat.trial   = singleTrials;
  stat.dimord  = 'trial_chan_chan_time';
end
stat.label       = spike.label(chansel);

% do the general cleanup and bookkeeping at the end of the function

ft_postamble previous   spike
ft_postamble provenance stat
ft_postamble history    stat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C] = spike_crossx_matlab(tX,tY,binsize,nbins)
tX = sort(tX(:));
tY = sort(tY(:));

minLag = - binsize * (nbins-1) / 2;
j = 0:nbins-1;
B = minLag + j * binsize;
tX(tX<(tY(1)+minLag) | tX>(tY(end)-minLag))   = [];
if isempty(tX)
  C = zeros(1,length(B)-1);
  return;
end
tY(tY>(tX(end)-minLag) | tY<(tX(1)+minLag)) = [];
if isempty(tY)
  C = zeros(1,length(B)-1);
  return;
end
nX = length(tX); nY = length(tY);

% compute all distances at once using a multiplication trick
if (nX*nY)<2*10^7  % allow matrix to grow to about 150 MB, should always work
  D = log(exp(-tX(:))*exp(tY(:)'));
  D = D(:);
  D(abs(D)>abs(minLag)) = [];
  [C] = histc(D,B);
  C(end) = [];
else
  % break it down in pieces such that nX*nY<2*10*7
  k   = 2;
  nXs = round(nX/k);
  while (nXs*nY)>=(2*10^7)
    k = k+1;
    nXs = round(nX/k);
  end
  
  % get the indices
  steps = round(nX/k);
  begs = 1:steps:steps*k;
  ends = begs+(steps-1);
  rm   = begs>nX;
  begs(rm) = [];
  ends(rm) = [];
  ends(end) = nX;
  nSteps    = length(begs);
  
  C = zeros(1,length(B));
  for iStep = 1:nSteps
    d = log(exp(-tX(begs(iStep):ends(iStep)))*exp(tY(:)'));
    d = d(:);
    d(abs(d)>abs(minLag)) = [];
    C = C + histc(d,B)';
  end
  C(end) = [];
end

if isempty(C)
  C = zeros(1,length(B)-1);
end
