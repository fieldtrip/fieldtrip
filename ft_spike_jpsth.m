function [jpsth] = jointpsth(cfg,psth)

% FT_spike_JPSTH computes the joint peristimulus histograms for spiketrains and a shift
% predictor (e.g. see Aertsen et al. 1989).
%
% Use as
%   [JPSTH] = FT_spike_JPSTH(CFG,PSTH)
%
% The shift predictor is computed in consecutive trials in a symmetric way. For
% example, we compute the jpsth for chan 1 in trial 1 versus chan 2 in trial 2, but also
% for chan 1 in trial 2 versus chan 2 in trial 1. This gives (nTrials-1)*2 jpsth matrices
% for individual trials. Picking consecutive trials and computing the shift predictor in 
% a symmetric way ensures that slow changes in the temporal structure do not affect the shift 
% predictor. 
%
% The input PSTH should be organised as the input from FT_SPIKE_PSTH, FT_SPIKE_DENSITY or
% TIMELOCKANALYSIS containing a field PSTH.trial and PSTH.time 
% In any case, one is expected to use cfg.keeptrials = 'yes' in these functions
%   
% Configurations options:
%
%   cfg.normalization    = 'no' (default), or 'yes'.  If requested (see cfg.normalization), the joint 
%                           psth is normalized as in van Aertsen et al. (1989), by 
%                           D(u,v)/sqrt(D(u,u)*D(v,v) w
%                           here D(u,v) is the difference of the average jpsth with the predicted jpsth 
%                           (see ref. for details), giving a quantity between -1 and 1. 
%                           Since this method normalizes by the mean across all trials, it can be
%                           confounded by latency drifs over trials.                         
%   cfg.shiftpredictor   = 'no' (default) or 'yes'. If 'yes', then JPSTH.AVG, JPSTH.VAR
%                           AND JPSTH.DOF will apply to the shiftpredictor, not the jpsth
%                           proper.                           
%   cfg.channelcmb       =  Mx2 cell-array with selection of channel pairs (default = {'all' 'all'}),
%                           see FT_CHANNELCOMBINATION for details
%   cfg.trials           = 'all' (default) or numerical or logical array of to be selected trials.
%   cfg.latency          = [begin end] in seconds, 'maxperiod' (default), 'prestim'(t<=0), or
%                          'poststim' (t>=0).                          
%   cfg.keeptrials       = 'yes' or 'no' (default).
%

% Copyright (C) 2010, Martin Vinck

% check the configuration inputs and enter the defaults
defaults.shiftpredictor = {'no' 'yes'};          
defaults.normalization  = {'no' 'yes'};            
defaults.trials         = {'all'};              
defaults.latency        = {'maxperiod'};              
defaults.keeptrials     = {'no' 'yes'};               
defaults.channelcmb     = {{'all','all'}};              
cfg = ft_spike_sub_defaultcfg(cfg,defaults);

% check whether PSTH contains single trials
if ~all(isfield(psth,{'trial' 'time'})), 
  error('MATLAB:ft_spike_jpsth:psth:missingFields',...
        'input PSTH should contain the field trial and time')
end

% get the number of trials or change DATA according to cfg.trials
if  strcmp(cfg.trials,'all')    
    cfg.trials = 1:size(psth.trial,1);
elseif islogical(cfg.trials)
    cfg.trials = find(cfg.trials); 
elseif ~isrealvec(cfg.trials);
    error('MATLAB:ft_spike_jpsth:cfg:trials:wrongInput',...
    'cfg.trials should be logical or numerical selection or string "all"');  
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
   if maxTime<=0, error('MATLAB:ft_spike_jpsth:cfg:latency',...
    'cfg.latency = "poststim" only allowed if psth.time(end)>0')
   end   
elseif strcmp(cfg.latency,'prestim')
   if minTime>=0, error('MATLAB:ft_spike_jpsth:cfg:latency',...
    'cfg.latency = "prestim" only allowed if psth.time(1)<0')
   end
   cfg.latency = [minTime 0]; %seems fishy, what if minTime > 0? CHECK OTHER FUNCS AS WELL
elseif ~isrealvec(cfg.latency)||length(cfg.latency)~=2
  error('MATLAB:ft_spike_jpsth:cfg:latency',...
    'cfg.latency should be "max" or 1-by-2 numerical vector');
end
if cfg.latency(1)>=cfg.latency(2),
     error('MATLAB:ft_spike_jpsth:cfg:latency:wrongOrder',...
    'cfg.latency(2) should be greater than cfg.latency(1)')
end
% check whether the time window fits with the data
if (cfg.latency(1) < minTime), cfg.latency(1) = minTime; 
  warning('MATLAB:ft_spike_jpsth:correctLatencyBeg',...
          'Correcting begin latency of averaging window');
end
if (cfg.latency(2) > maxTime), cfg.latency(2) = maxTime;
  warning('MATLAB:ft_spike_jpsth:correctLatencyEnd',...
          'Correcting end latency of averaging window');
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
if nCmbs==0, error('MATLAB:ft_spike_jpsth:cfg:channelcmb:noneSelected', ...
    'No channel combination selected')
end
  
% decompose into single channels 
chanSel = unique(cmbindx(:)); % this gets sorted ascending by default
nChans  = length(chanSel);

% preallocate avg in chan x chan format, this can take more memory, but its more intuitive
if strcmp(cfg.keeptrials,'yes'), singleTrials = zeros(nTrials,nBins,nBins,nCmbs); end
[avgJpsth,varJpsth] = deal(zeros(nBins,nBins,nChans,nChans)); 
if strcmp(cfg.shiftpredictor,'yes'), 
  [avgshiftpredictor,varshiftpredictor] = deal(zeros(nBins,nBins,nChans,nChans)); 
end

% compute the dof for the jpsth, recompute this since we selected trials.
sumNum = sum(~isnan(psth.trial(:,1,:)));
dof    = squeeze(sumNum)';
dof1   = dof(ones(nBins,1),:);
dof    = min(dof1,dof1'); % the actual dof is always given by the maximum of the 2

% check if we have to compute dof while in the loop 
computeDofShift = any(dof(:)~=size(psth.trial,1)); %i.e, we have nans in psth?
if computeDofShift
  dofShift = zeros(nBins); % will compute this on the fly
else
  dofShift = 2*(nTrials-1);
end

% compute the joint psth
for iCmb = 1:nCmbs
    indx1 = cmbindx(iCmb,1);
    indx2 = cmbindx(iCmb,2);
    [ss,s]  = deal(sparse(nBins,nBins));

    if strcmp(cfg.shiftpredictor,'yes'), [ssShift,sShift]  = deal(s); end

    % already compute the quantities to normalize the jpsth
    if strcmp(cfg.normalization,'yes')
      meanX = squeeze(nanmean(psth.trial(:,indx1,:))); % psth can contain nans
      repMeanX = meanX(:,ones(1,nBins));    
      meanH = squeeze(nanmean(psth.trial(:,indx2,:)))';
      repMeanH = meanH(ones(nBins,1),:);
      meanXH   = repMeanX.*repMeanH;

      % compute the product of variances
      varX = squeeze(nanvar(psth.trial(:,indx1,:),1,1));
      varH = squeeze(nanvar(psth.trial(:,indx2,:),1,1))';
      
      % var is numerically instable when var=0, check if diff is different from 0
      % this way we can set the jpsth to correct value of 0 later
      d = diff(psth.trial(:,indx1,:),[],1);
      hasDiff = squeeze(any(d,1));
      varX(~(hasDiff>0)) = 0;      
      d = diff(psth.trial(:,indx2,:),[],1);
      hasDiff = squeeze(any(d,1));
      varH(~(hasDiff>0)) = 0;

      repVarX = varX(:,ones(1,nBins));
      repVarH = varH(ones(nBins,1),:);
      varXH   = repVarX.*repVarH;
            
    end

    for iTrial = 1:nTrials
        
        x = squeeze(psth.trial(iTrial,indx1, :)); % first chan
        h = squeeze(psth.trial(iTrial,indx2, :)); % second chan
        
        % get the indices of the channels
        indX = find(x>0).';
        indH = find(h>0); 
        nX = length(indX);
        nH = length(indH);

        if computeDofShift % otherwise we do not have nans
          nanBins         = isnan(x)';
          currentNan      = nanBins(ones(nBins,1),:);
          currentTrialNan = currentNan | currentNan';
        end

        % compute the shift predictor symmetrically, x21 to x12, and x22 to x11
        if strcmp(cfg.shiftpredictor,'yes') && iTrial>1
          
          % in case of shift predictor, compute the dof, is same for all cmbs
          if computeDofShift && iCmb==1 % otherwise we do not have nans
            nanBinsOld = nanBinsOld';
            oldNan     = nanBinsOld(:,ones(1,nBins));
            hasNum     = ~currentNan & ~oldNan;         
            dofShift   = dofShift + double(hasNum + hasNum.');
          end
          
          % might occur in neighbouring trials
          nX_old = length(indX_old);
          nH_old = length(indH_old);

          % first from x21 to x12
          repX = indX(ones(nH_old,1),:);
          repX = repX(:);
          repH = indH_old(:,ones(1,nX));
          repH = repH(:);        
          vals = x(repX).*hOld(repH); % calculate the values for every bin

          % construct the psth matrix per trial and add to the sum
          jpsthTrial = sparse(repX,repH,vals,length(x),length(h));
          sShift  = sShift  + jpsthTrial; % now dof goes times 2
          ssShift = ssShift + jpsthTrial.^2;
        
          % then from x22 to x11
          repX = indX_old(ones(nH,1),:);
          repX = repX(:);
          repH = indH(:,ones(1,nX_old));
          repH = repH(:);                            
          vals = xOld(repX).*h(repH); % calculate the values for every bin
 
          % construct the psth matrix per trial, 
          jpsthTrial = sparse(repX,repH,vals,length(x),length(h));        
          sShift = sShift   + jpsthTrial; % now dof goes times 2
          ssShift = ssShift + jpsthTrial.^2;
        end

        if strcmp(cfg.shiftpredictor,'yes')
          indH_old = indH;
          indX_old = indX;
          xOld = x;
          hOld = h;
          if computeDofShift,  nanBinsOld = nanBins; end
        end
        
        % replicate the indices matrices into two grids and vectorize again
        repH = indH(:,ones(1,nX));
        repH = repH(:);        
        repX = indX(ones(nH,1),:);
        repX = repX(:);
        vals = x(repX).*h(repH); % calculate the values for every bin
        jpsthTrial = sparse(repX,repH,vals,length(x),length(h));  % construct the psth matrix per trial, 
        
        % compute the sum and the squared sum (for the variance)
        s  = s  + jpsthTrial;
        ss = ss + jpsthTrial.^2;                            
        if strcmp(cfg.keeptrials,'yes'), 
           if computeDofShift, jpsthTrial(currentTrialNan) = NaN; end % store as Nans in single trials.
           singleTrials(iTrial,:,:,iCmb) = jpsthTrial; 
        end
    end

     % get the individual channel numbers (not in original data.label sense)
    indx1 = find(chanSel==indx1);
    indx2 = find(chanSel==indx2);

    % compute the raw average and the raw variance
    m = s./dof; % still delete the 0 dof
    if strcmp(cfg.normalization,'yes')
      m = (m - meanXH) ./ sqrt(varXH); 
      m(meanXH==0) = 0; % with no spikes in joint bin there, jpsth should be 0        
      m(varXH==0)  = 0;
    end
    m(dof==0) = NaN; % should we output nans were there was no data? what for psth?
    
    avgJpsth(:,:,indx1,indx2)       = m;
    avgJpsth(:,:,indx2,indx1)       = m';
    
    v = (ss - s.^2./dof)./(dof-1);
    v(dof==0) = NaN;
    varJpsth(:,:,indx1,indx2)       = v;
    varJpsth(:,:,indx2,indx1)       = v';

    % compute the average and the variance for the shift predictor.
    if strcmp(cfg.shiftpredictor,'yes')
      m = sShift ./ dofShift;
      if strcmp(cfg.normalization,'yes')
        m = (m - meanXH) ./ sqrt(varXH); 
        m(meanXH==0) = 0; % with no spikes in joint bin there, jpsth should be 0
        m(varXH==0)  = 0;
      end
      m(dof==0) = NaN;
      avgshiftpredictor(:,:,indx1,indx2)       = m;
      avgshiftpredictor(:,:,indx2,indx1)       = m';    
      
      v = (ssShift - sShift.^2./dofShift)./(dofShift-1);
      v(dof==0) = NaN;
      varshiftpredictor(:,:,indx1,indx2)       = v;
      varshiftpredictor(:,:,indx2,indx1)       = v';
    end         
end
   
% collect the results
if strcmp(cfg.shiftpredictor,'no')
  jpsth.avg      = avgJpsth;
  jpsth.var      = varJpsth;
  jpsth.dof      = dof;
else
  jpsth.avg       = avgshiftpredictor;
  jpsth.var      = varshiftpredictor;
  jpsth.dof      = dofShift;  
end
jpsth.time     = psth.time;
try jpsth.bins = psth.bins; end
jpsth.label    = psth.label(chanSel); % keep this as reference for JPSTH.avg
jpsth.labelcmb = psth.label(cmbindx); % keep this as reference which channel cmbs have values
jpsth.psth     = psth.avg(chanSel,:);
if (strcmp(cfg.keeptrials,'yes'))
  jpsth.trial = singleTrials;
  jpsth.dimord = 'rpt_time_time_chan';
else
  jpsth.dimord = 'time_time_chan_chan';
end

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
% remember the configuration details of the input data
if isfield(psth,'cfg'), cfg.previous = psth.cfg; end
% remember the exact configuration details in the output 
jpsth.cfg     = cfg;


