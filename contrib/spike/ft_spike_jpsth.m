function [stat] = ft_spike_jpsth(cfg, psth)

% FT_SPIKE_JPSTH computes the joint peristimulus histograms for spiketrains
% and a shift predictor (for example see Aertsen et al. 1989).
%
% The shift predictor is computed in consecutive trials in a symmetric way.
% For example, we compute the jpsth for chan 1 in trial 1 versus chan 2 in
% trial 2, but also for chan 1 in trial 2 versus chan 2 in trial 1. This
% gives (nTrials-1)*2 jpsth matrices for individual trials. Picking
% consecutive trials and computing the shift predictor in a symmetric way
% ensures that slow changes in the temporal structure do not affect the
% shift predictor (as opposed to shuffling the order of all trials for one
% of the two channels).
%
% Use as
%   [jpsth] = ft_spike_jpsth(cfg,psth)
%
% The input PSTH should be organised as the input from FT_SPIKE_PSTH,
% FT_SPIKE_DENSITY or FT_TIMELOCKANALYSIS containing a field PSTH.trial and
% PSTH.time. In any case, one is expected to use cfg.keeptrials = 'yes' in
% these functions.
%
% Configurations:
%   cfg.method           = 'jpsth' or 'shiftpredictor'. If 'jpsth', we
%                           output the normal stat. If 'shiftpredictor',
%                           we compute the jpsth after shuffling subsequent
%                           trials.
%   cfg.normalization    = 'no' (default), or 'yes'. If requested, the joint psth is normalized as in van Aertsen et al. (1989).
%   cfg.channelcmb       =  Mx2 cell-array with selection of channel pairs (default = {'all' 'all'}), see FT_CHANNELCOMBINATION for details
%   cfg.trials           = 'all' (default) or numerical or logical array of to be selected trials.
%   cfg.latency          = [begin end] in seconds, 'maxperiod' (default), 'prestim'(t<=0), or 'poststim' (t>=0)
%   cfg.keeptrials       = 'yes' or 'no' (default)
%
% See also FT_SPIKE_PSTH

% Copyright (C) 2010, Martin Vinck
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
ft_preamble provenance psth
ft_preamble trackconfig

psth = ft_checkdata(psth, 'datatype', 'timelock', 'hastrials', 'yes', 'feedback', 'yes');

% get the default options
cfg.trials         = ft_getopt(cfg,'trials', 'all');
cfg.latency        = ft_getopt(cfg,'latency','maxperiod');
cfg.keeptrials     = ft_getopt(cfg,'keeptrials', 'no');
cfg.method         = ft_getopt(cfg,'method', 'jpsth');
cfg.normalization  = ft_getopt(cfg,'normalization', 'no');
cfg.channelcmb     = ft_getopt(cfg,'channelcmb', 'all');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'latency', {'char', 'ascendingdoublebivector'});
cfg = ft_checkopt(cfg,'trials', {'char', 'doublevector', 'logical'});
cfg = ft_checkopt(cfg,'keeptrials', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'method', 'char', {'jpsth', 'shiftpredictor'});
cfg = ft_checkopt(cfg,'normalization', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'channelcmb', {'char', 'cell'});

% reject configuration inputs that are not processed
cfg = ft_checkconfig(cfg, 'allowed', {'latency', 'trials', 'keeptrials', 'method', 'normalization', 'channelcmb'});

% get the number of trials or change DATA according to cfg.trials
if  strcmp(cfg.trials,'all')
  cfg.trials = 1:size(psth.trial,1);
elseif islogical(cfg.trials) || all(cfg.trials==0 | cfg.trials==1)
  cfg.trials = find(cfg.trials);
end
cfg.trials  = sort(cfg.trials(:));
psth.trial  = psth.trial(cfg.trials,:,:);
nTrials     = length(cfg.trials);

% select the time
minTime     = psth.time(1);
maxTime     = psth.time(end); % we know this is ordered vector
if strcmp(cfg.latency,'maxperiod')
  cfg.latency = [minTime maxTime];
elseif strcmp(cfg.latency,'poststim')
  cfg.latency = [0 maxTime];
  if maxTime<=0, error('cfg.latency = "poststim" only allowed if psth.time(end)>0'); end
elseif strcmp(cfg.latency,'prestim')
  if minTime>=0, error('cfg.latency = "prestim" only allowed if psth.time(1)<0'); end
  cfg.latency = [minTime 0]; %seems fishy, what if minTime > 0? CHECK OTHER FUNCS AS WELL
end
% check whether the time window fits with the data
if (cfg.latency(1) < minTime), cfg.latency(1) = minTime;
  warning('Correcting begin latency of averaging window');
end
if (cfg.latency(2) > maxTime), cfg.latency(2) = maxTime;
  warning('Correcting end latency of averaging window');
end

% get the right indices in psth.time and select this part of the data
indx       = nearest(psth.time,cfg.latency(1)) : nearest(psth.time,cfg.latency(2));
psth.time  = psth.time(indx);
psth.trial = psth.trial(:,:,indx);
nBins      = length(psth.time);

% determine the corresponding indices of the requested channel combinations
cfg.channelcmb = ft_channelcombination(cfg.channelcmb, psth.label);
cmbindx = zeros(size(cfg.channelcmb));
for k=1:size(cfg.channelcmb,1)
  cmbindx(k,1) = strmatch(cfg.channelcmb(k,1), psth.label, 'exact');
  cmbindx(k,2) = strmatch(cfg.channelcmb(k,2), psth.label, 'exact');
end
nCmbs 	   = size(cmbindx,1);
if nCmbs==0, error('No channel combination selected'); end

% decompose into single channels
chanSel = unique(cmbindx(:)); % this gets sorted ascending by default
nChans  = length(chanSel);

% preallocate avg in chan x chan format, this can take more memory, but its more intuitive
if strcmp(cfg.keeptrials,'yes')
  singleTrials = NaN(nTrials,nChans,nChans,nBins,nBins);
  warning('storing single trials for jpsth is memory expensive, please check');
end
[out,varOut,dofOut] = deal(zeros(nChans,nChans,nBins,nBins));

% compute the joint psth
ft_progress('init', 'text', 'Please wait...');
for iCmb = 1:nCmbs
  indxData1 = cmbindx(iCmb,1); % index for the data
  indxData2 = cmbindx(iCmb,2);
  indxOut1 = find(chanSel==indxData1); % this is in the order of the output
  indxOut2 = find(chanSel==indxData2);
  
  [ss,s,df]  = deal(zeros(nBins,nBins));
  
  % already compute the quantities to normalize the jpsth
  if strcmp(cfg.normalization,'yes')
    mean1    = squeeze(nanmean(psth.trial(:,indxData1,:))); % psth can contain nans
    mean2    = squeeze(nanmean(psth.trial(:,indxData2,:)))';
    mean12   = mean1(:)*mean2(:)';
    diff1    = nansum(diff(squeeze(psth.trial(:,indxData1,:)),[],1),1); % this is just to avoid rounding errors, as var gives these
    diff2    = nansum(diff(squeeze(psth.trial(:,indxData2,:)),[],1),1); % this is just to avoid rounding errors, as var gives these
    var1     = squeeze(nanvar(psth.trial(:,indxData1,:),1,1));
    var1(diff1==0) = 0;
    var2     = squeeze(nanvar(psth.trial(:,indxData2,:),1,1))';
    var2(diff2==0) = 0;
    var12    = var1(:)*var2(:)';
    var12(mean12==0) = 0;
  end
  
  for iTrial = 1:nTrials
    ft_progress(iTrial/nTrials, 'Processing trial %d from %d for combination %d out of %d', iTrial, nTrials, iCmb, nCmbs);
    psth1 = squeeze(psth.trial(iTrial,indxData1, :)); % first chan
    psth2 = squeeze(psth.trial(iTrial,indxData2, :)); % second chan
    isNum1 = double(~isnan(psth1));
    isNum2 = double(~isnan(psth2));
    
    switch cfg.method
      case 'jpsth'
        
        % compute the 2-D product with the matrix multiplication
        jpsthTrial = psth1(:)*psth2(:)';
        dofTrial   = isNum1(:)*isNum2(:)';
        
        % compute the sum and the squared sum (for the variance)
        s(dofTrial>0)  = s(dofTrial>0)  + jpsthTrial(dofTrial>0);
        ss(dofTrial>0) = ss(dofTrial>0) + jpsthTrial(dofTrial>0).^2;
        df(dofTrial>0) = df(dofTrial>0) + dofTrial(dofTrial>0);
        jpsthTrial(dofTrial==0) = NaN;
        if strcmp(cfg.keeptrials,'yes')
          singleTrials(iTrial,indxOut1,indxOut2,:,:) = jpsthTrial;
          singleTrials(iTrial,indxOut2,indxOut1,:,:) = jpsthTrial';
        end
        
      case 'shiftpredictor'
        
        if iTrial>1
          psth1Prev    = squeeze(psth.trial(iTrial-1,indxData1, :)); % first chan
          psth2Prev    = squeeze(psth.trial(iTrial-1,indxData2, :)); % second chan
          isNum1Prev = double(~isnan(psth1Prev));
          isNum2Prev = double(~isnan(psth2Prev));
          
          jpsthTrial = nansum(cat(3,psth1(:)*psth2Prev(:)',psth1Prev(:)*psth2(:)'),3);
          dofTrial   = nansum(cat(3,isNum1(:)*isNum2Prev(:)',isNum1Prev(:)*isNum2(:)'),3);
          
          s(dofTrial>0)     = s(dofTrial>0)   + jpsthTrial(dofTrial>0); % now dof goes times 2
          ss(dofTrial>0)    = ss(dofTrial>0)  + jpsthTrial(dofTrial>0).^2;
          df(dofTrial>0)    = df(dofTrial>0) + dofTrial(dofTrial>0);
          jpsthTrial(dofTrial==0) = NaN;
          jpsthTrial = jpsthTrial./dofTrial; % normalize for having two combinations
          
          if strcmp(cfg.keeptrials,'yes')
            singleTrials(iTrial,indxOut1,indxOut2,:,:) = jpsthTrial;
            singleTrials(iTrial,indxOut2,indxOut1,:,:) = jpsthTrial';
          end
        end
    end
  end
  
  % compute the mean and the variance of the output
  m = s./df; % still delete the 0 dof
  if strcmp(cfg.normalization,'yes')
    m = (m - mean12) ./ sqrt(var12);
    m(mean12==0) = 0; % with no spikes in joint bin there, jpsth should be 0
    m(var12==0)  = 0; % if variance is zero, we assume 0/0 = 0
  end
  m(df==0) = NaN; % no trials: must be a NaN
  
  out(indxOut1,indxOut2,:,:)       = m;
  out(indxOut2,indxOut1,:,:)       = m';
  
  v = (ss - s.^2./df)./(df-1);
  v(df<=1) = NaN;
  varOut(indxOut1,indxOut2,:,:)    = v;
  varOut(indxOut2,indxOut1,:,:)    = v';
  
  dofOut(indxOut1,indxOut2,:,:)    = df;
  dofOut(indxOut2,indxOut1,:,:)    = df';
  
end % for iCmb
ft_progress('close')

% collect the results
if strcmp(cfg.method,'jpsth')
  stat.jpsth      = out;
else
  stat.shiftpredictor = out;
end
stat.var        = varOut;
stat.dof        = df;
stat.time       = psth.time;
stat.psth       = shiftdim(mean(psth.trial(:,chanSel,:), 1), 1); % the input is single-trials, compute the mean over selected trials
stat.label      = psth.label(chanSel); % keep this as reference for JPSTH.avg
if (strcmp(cfg.keeptrials,'yes'))
  stat.trial = singleTrials;
  stat.dimord = 'rpt_time_time_chan_chan';
else
  stat.dimord = 'time_time_chan_chan';
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble previous   psth
ft_postamble provenance stat
ft_postamble history    stat
