function [stat] = ft_spike_phaselockstat(cfg,sts)

% FT_SPIKE_PHASELOCKSTAT computes phase-locking statistics for spike-LFP
% phases. These contain the PPC statistics published in Vinck et al. 2010
% (Neuroimage) and Vinck et al. 2011 (Journal of Computational
% Neuroscience).
%
% Use as:
%   [STAT = FT_SPIKE_PHASELOCKSTAT(CFG,STS)
%
% Inputs:
%   STS should be a structure as obtained from from the FT_SPIKE_TRIGGEREDSPECTRUM or FT_SPIKETRIGGEREDSPECTRUM function. 
%
% Configurations (CFG) related to calculation of phaselockingvalue
%
% cfg.method  = string, indicating which statistic to compute. Can be:
%     -'plv' : phase-locking value, computes the resultant length over spike
%              phases. Is positively biased by sample size.
%     -'ang' : computes the angular mean of the spike phases.
%     -'ral' : computes the rayleigh statistic.
%     -'ppc0': computes the pairwise-phase consistency across all available
%              spike pairs (Vinck et al., 2010).
%     -'ppc1': computes the pairwise-phase consistency across all available
%              spike pairs with exclusion of spike pairs in the same trial.
%              This avoids history effects within spike-trains to influence
%              phase lock statistics.
%     -'ppc2': computes the PPC across all spike pairs with exclusion of
%              spike pairs in the same trial, but applies a normalization
%              for every set of trials. This estimator has more variance but
%              is more robust against dependencies between spike phase and
%              spike count.
%         
% cfg.timwin  = double or 'all' (default)
%   - doube: indicates we compute statistic with a
%            sliding window of cfg.timwin, i.e. time-resolved analysis.
%   - 'all': we compute statistic over all time-points,
%            i.e. in non-time resolved fashion.
%
% cfg.winstepsize  = double, stepsize of sliding window in seconds. For
%   example if cfg.winstepsize = 0.1, we compute stat every other 100 ms.
%
% cfg.channel      = Nx1 cell-array or numerical array with selection of
%   channels (default = 'all'),See CHANNELSELECTION for details
%
% cfg.spikesel     = 'all' (default) or numerical or logical selection of spikes.
%
% cfg.foi          = 'all' or numerical vector that specifies a subset of
%   frequencies in Hz (and not in indices!).                                    
%
% cfg.latency      = [beg end] in sec, or 'maxperiod',  'poststim' or
%  'prestim'.  This determines the start and end of time-vector.
%
% cfg.spikechannel = label of ONE unit, according to FT_CHANNELSELECTION
%
% cfg.chanavg      = 'weighted', 'unweighted' (default) or 'no'.
%  - 'weighted'  : we average across channels by weighting by the LFP power.
%                  This is identical to adding the raw LFP signals in time 
%                  and then taking their FFT.
%  - 'unweighted': we average across channels after normalizing for LFP power. 
%                  This is identical to filtering LFP signals, normalizing 
%                  their power, averaging them, and then taking their FFT.
%  - 'no'        : no weighting is performed, statistic is computed for
%                  every LFP channel.
% cfg.trials       = vector of indices (e.g., 1:2:10),
%                    logical selection of trials (e.g., [1010101010]), or
%                   'all' (default)
%
% Main outputs:
%     stat.doftrial                   =  nTimepoints-by-nChan-by-nFreqs number of trials
%     stat.dofspike                   =  nTimepoints-by-nChan-by-nFreqs number of spikes
%     stat.label                      =  name of unit;
%     stat.(cfg.method)               =  nChan-by-nFreqs statistic
%     stat.lfplabel                   =  nChans cell array with LFP labels
%

%   Copyright (c) Martin Vinck (2012)
%
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% check if the data is of sts format, and convert from old format if required
sts = ft_checkdata(sts,'datatype', 'spike', 'feedback', 'yes');

% get the options
cfg.method         = ft_getopt(cfg,'method', 'ppc1');
cfg.channel        = ft_getopt(cfg,'channel', 'all');
cfg.spikechannel   = ft_getopt(cfg,'spikechannel', sts.label{1});
cfg.latency        = ft_getopt(cfg,'latency', 'maxperiod');
cfg.spikesel       = ft_getopt(cfg,'spikesel', 'all');
cfg.chanavg        = ft_getopt(cfg,'chanavg', 'no');
cfg.foi            = ft_getopt(cfg,'foi', 'all');
cfg.trials         = ft_getopt(cfg,'trials', 'all');
cfg.timwin         = ft_getopt(cfg,'timwin', 'all');
cfg.winstepsize    = ft_getopt(cfg,'winstepsize', 0.01);

% ensure that the options are valid
cfg = ft_checkopt(cfg,'method', 'char', {'ppc0', 'ppc1', 'ppc2', 'ang', 'ral', 'plv'});
cfg = ft_checkopt(cfg,'foi',{'char', 'double'});
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'channel', {'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'spikesel', {'char', 'logical', 'double'});
cfg = ft_checkopt(cfg,'chanavg', 'char', {'weighted', 'unweighted', 'no'});
cfg = ft_checkopt(cfg,'latency', {'char', 'doublevector'});
cfg = ft_checkopt(cfg,'trials', {'char', 'double', 'logical'}); 
cfg = ft_checkopt(cfg,'timwin', {'double', 'char'}); 
cfg = ft_checkopt(cfg,'winstepsize', {'double'}); 

cfg = ft_checkconfig(cfg,'allowed', {'method', 'channel', 'spikechannel', 'latency', 'spikesel', 'chanavg', 'foi', 'trials', 'timwin', 'winstepsize'});

% collect channel information
cfg.channel        = ft_channelselection(cfg.channel, sts.lfplabel);
chansel            = match_str(sts.lfplabel, cfg.channel); 

% get the spikechannels
spikelabel       = sts.label;
cfg.spikechannel = ft_channelselection(cfg.spikechannel, spikelabel);
unitsel          = match_str(spikelabel, cfg.spikechannel);
nspikesel        = length(unitsel); % number of spike channels
if nspikesel>1, error('only one unit should be selected for now'); end

% collect frequency information
if strcmp(cfg.foi, 'all'),  
  cfg.foi           = sts.freq; 
  freqindx          = 1:length(sts.freq);
else
  freqindx  = zeros(1,length(foi));
  for iFreq = 1:length(cfg.foi)
    freqindx(iFreq)         = nearest(sts.freq,cfg.foi(iFreq)); 
  end
end
if length(freqindx)~=length(unique(freqindx)) 
  error('Please select every frequency only once, are you sure you selected in Hz?')
end
nFreqs     = length(freqindx);
cfg.foi    = sts.freq(freqindx); % update the information again

% create the spike selection 
nSpikes = length(sts.trial{unitsel});
if strcmp(cfg.spikesel,'all'), 
  cfg.spikesel = 1:length(sts.trial{unitsel});
elseif islogical(cfg.spikesel)
  cfg.spikesel = find(cfg.spikesel);
end
if ~isempty(cfg.spikesel)
  if max(cfg.spikesel)>nSpikes || length(cfg.spikesel)>nSpikes 
    error('cfg.spikesel must not exceed number of spikes and select every spike just once')
  end
end
  
% select on basis of latency
if strcmp(cfg.latency, 'maxperiod'),
   cfg.latency = [min(sts.trialtime(:)) max(sts.trialtime(:))];
elseif strcmp(cfg.latency,'prestim')
   cfg.latency = [min(sts.trialtime(:)) 0];
elseif strcmp(cfg.latency,'poststim')
   cfg.latency = [0 max(sts.trialtime(:))];
elseif ~isrealvec(cfg.latency) && length(cfg.latency)~=2 
   error('cfg.latency should be "maxperiod", "prestim", "poststim" or 1-by-2 numerical vector');
end
if cfg.latency(1)>=cfg.latency(2), 
   error('cfg.latency should be a vector in ascending order, i.e., cfg.latency(2)>cfg.latency(1)');
end
inWindow = find(sts.time{unitsel}>=cfg.latency(1) & cfg.latency(2)>=sts.time{unitsel});

% selection of the trials
cfg        = trialselection(cfg,sts);

% do the final selection, and select on sts structure
isintrial    = ismember(sts.trial{unitsel}, cfg.trials);
spikesel     = intersect(cfg.spikesel(:),inWindow(:));
spikesel     = intersect(spikesel,find(isintrial));
cfg.spikesel = spikesel; %intersect(spikesel(:),isNum(:));
spikenum     = length(cfg.spikesel); % number of spikes that were finally selected
if isempty(spikenum), warning('No spikes were selected after applying cfg.latency, cfg.spikesel and cfg.trials'); end
sts.fourierspctrm = sts.fourierspctrm{unitsel}(cfg.spikesel,chansel,freqindx);
sts.time          = sts.time{unitsel}(cfg.spikesel);
sts.trial         = sts.trial{unitsel}(cfg.spikesel);

% average the lfp channels (weighted, unweighted, or not)
if strcmp(cfg.chanavg,'unweighted')
  sts.fourierspctrm = sts.fourierspctrm ./ abs(sts.fourierspctrm); % normalize the angles before averaging   
  sts.fourierspctrm = nanmean(sts.fourierspctrm,2); % now rep x 1 x freq
  nChans = 1;
elseif strcmp(cfg.chanavg,'no')
  nChans = length(chansel);
elseif strcmp(cfg.chanavg,'weighted')
  sts.fourierspctrm = nanmean(sts.fourierspctrm,2); % now rep x 1 x freq
  nChans = 1;
end

% normalize the spectrum first
sts.fourierspctrm = sts.fourierspctrm ./ abs(sts.fourierspctrm); % normalize the angles before averaging   

if strcmp(cfg.timwin,'all')

  switch cfg.method
    case 'ang'
     out  = angularmean(sts.fourierspctrm);
    case 'plv'
     out  = resultantlength(sts.fourierspctrm);
    case 'ral'
     out  = rayleightest(sts.fourierspctrm); % the rayleigh test
    case 'ppc0'
     out  = ppc(sts.fourierspctrm);
    case {'ppc1', 'ppc2'}
      
      % check the final set of trials present in the spikes
      trials = unique(sts.trial);

      % loop init for PPC 2.0
      [S,SS,dof,dofSS] = deal(zeros(1,nChans,nFreqs));
      nTrials = length(trials);
      for iTrial = 1:nTrials % compute the firing rate
          trialNum      = trials(iTrial);
          spikesInTrial = find(sts.trial == trialNum);
          spcU          = sts.fourierspctrm(spikesInTrial,:,:);
          spc           = spcU./abs(spcU);

          % compute PPC 2.0 according to Vinck et al. (2011) using summation per
          % trial
          if strcmp(cfg.method,'ppc1')
            if ~isempty(spc)
              m = nanmean(spc,1); % no problem with NaN
              hasNum = ~isnan(m);
              S(hasNum)  = S(hasNum)  + m(hasNum); % add the imaginary and real components
              SS(hasNum) = SS(hasNum) + m(hasNum).*conj(m(hasNum));
              dof(hasNum) = dof(hasNum) + 1; % dof needs to be kept per frequency
            end
          elseif strcmp(cfg.method,'ppc2')
            % compute PPC 1.0
            if ~isempty(spc)
               n      = sum(~isnan(spc),1);
               m      = nansum(spc,1); 
               hasNum = ~isnan(m);    
               S(hasNum)     = S(hasNum)  + m(hasNum); % add the imaginary and real components
               SS(hasNum)    = SS(hasNum) + m(hasNum).*conj(m(hasNum));
               dof(hasNum)   = dof(hasNum)  + n(hasNum);
               dofSS(hasNum) = dofSS(hasNum) + n(hasNum).^2;
            end                              
          end
      end

      [out]   = deal(NaN(1,nChans,nFreqs));
      hasTrl = dof>1;
      if strcmp(cfg.method,'ppc1')
        out(hasTrl) = (S(hasTrl).*conj(S(hasTrl)) - SS(hasTrl))./(dof(hasTrl).^2 - dofSS(hasTrl));
      else
        out(hasTrl) = (S(hasTrl).*conj(S(hasTrl)) - SS(hasTrl))./(dof(hasTrl).*(dof(hasTrl)-1));
      end
  end
  nSpikes = sum(~isnan(sts.fourierspctrm));
else % compute time-resolved spectra of statistic
  
  % make the sampling axis for the window
  bins      = cfg.latency(1):cfg.winstepsize:cfg.latency(2);
  N         = length(bins)-1; % number of bins
  wintime   = 0:cfg.winstepsize:cfg.timwin;
  win       = ones(1,length(wintime));
  stat.time = (bins(2:end)+bins(1:end-1))/2;
  if ~mod(length(win),2), win = [win 1]; end   % make sure the number of samples is uneven.

  out      = NaN(N,nChans,nFreqs);
  nSpikes  = zeros(N,nChans,nFreqs);
  for iChan = 1:nChans
    for iFreq = 1:nFreqs
  
      spctra       = sts.fourierspctrm(:,iChan,iFreq); % values to accumulate
      tm           = sts.time;
      hasnan       = isnan(spctra);    
      [tm(hasnan), spctra(hasnan)]  = deal([]);

      % compute the degree of freedom per time bin and the index for every spike
      [dof, indx] = histc(tm, bins); % get the index per spike, and number per bin                          
      if isempty(dof), continue,end    
      toDel       = indx==length(bins) | indx==0; % delete those points that are equal to last output histc or don't fall in
      [spctra(toDel),indx(toDel)] = deal([]);
      dof(end)    = []; % the last bin is a single point in time, so we delete it

      % compute the sum of spikes per window at every time point
      dof    = dof(:); % force it to be row
      dof    = conv2(dof(:),win(:),'same'); % get the total number of spikes across all trials
      nSpikes(:,iChan,iFreq) = dof;      
      
      switch cfg.method 
        case {'ang', 'plv', 'ppc0', 'ral'}

          % first create a vector with the phases at the samples
          x    = accumarray(indx(:),spctra,[N 1]); % simply the sum of the complex numbers
          
          % then compute the sum of the spectra for every timepoint
          y    = conv2(x(:),win(:),'same');     
          
          % now compute the output statistic
          hasnum  = dof>1;
          hasnum0 = dof>0;          
          if strcmp(cfg.method,'plv')
            out(hasnum,iChan,iFreq) = abs(y(hasnum)./dof(hasnum));
          elseif strcmp(cfg.method,'ral')
            Z = abs(y).^2./dof;
            P = exp(-Z) .* (1 + (2*Z - Z.^2)./(4*dof) -(24.*Z - 123*(Z.^2) + 76*(Z.^3) - 9*(Z.^4))./(288*dof.^2)); %Mardia 1972
            out(hasnum,iChan,iFreq) = P(hasnum);     
          elseif strcmp(cfg.method,'ang')
            out(hasnum0,iChan,iFreq) = angle(y(hasnum0)); %            
          elseif strcmp(cfg.method,'ppc0')
            out(hasnum,iChan,iFreq) = (y(hasnum).*conj(y(hasnum)) - dof(hasnum))./(dof(hasnum).*(dof(hasnum)-1)); % simplest form of ppc
          end

        case {'ppc1', 'ppc2'}
          trials  = unique(sts.trial);
          nTrials = length(trials);
          [S,SS,dofS,dofSS] = deal(zeros(length(bins)-1,1));
 
          % compute the new ppc versions
          for iTrial = 1:nTrials
            
             % select the spectra, time points, and trial numbers again
             trialNum      = trials(iTrial);
             spikesInTrial = find(sts.trial == trialNum);
             if isempty(spikesInTrial), continue,end
             spctraTrial  = sts.fourierspctrm(spikesInTrial,iChan,iFreq);       
             tm           = sts.time;
             hasnan       = isnan(spctraTrial);    
             [tm(hasnan), spctraTrial(hasnan)]  = deal([]);

             % bin the spikes and delete spikes out of the selected time
             [dof, indx]     = histc(tm(spikesInTrial), bins); % get the index per spike, and number per bin                          
             dof(end)        = [];
             toDel           = indx==length(bins) | indx==0; % delete those points that are equal to last output histc
             spctraTrial(toDel,:,:) = []; % delete those spikes from fourierspctrm as well
             indx(toDel) = []; % make sure index doesn't contain them
                          
             % first create a vector that sums the spctra
             x    = accumarray(indx(:),spctraTrial,[N 1]); 

             % then compute the moving average sum of this vector
             y    = conv2(x(:),win(:),'same'); % convolution again just means a sum                                  
             d    = conv2(dof(:),win(:),'same'); % get the dof of spikes per trial

             if strcmp(cfg.method,'ppc1')
               S     = S     + y; % add the imaginary and real components
               SS    = SS    + y.*conj(y);
               dofS  = dofS  + d(:);
               dofSS = dofSS + d(:).^2;
             else
               sl             = d(:)>0;
               m              = y./d(:);
               S(sl)          = S(sl)  + m(sl); % add the imaginary and real components
               SS(sl)         = SS(sl) + m(sl).*conj(m(sl));
               dofS(sl)       = dofS(sl) + 1;                     
             end
          end
          
          % we need at least two trials
          hasNum = dofS>1;
          if strcmp(cfg.method,'ppc1')
            out(hasNum,iChan,iFreq) = (S(hasNum).*conj(S(hasNum)) - SS(hasNum))./(dofS(hasNum).^2 - dofSS(hasNum));
          else
            out(hasNum,iChan,iFreq) = (S(hasNum).*conj(S(hasNum)) - SS(hasNum))./(dofS(hasNum).*(dofS(hasNum)-1));
          end
      end
    end
  end
end

% collect the outputs
outparam        = cfg.method;
stat.(outparam) = out;
stat.nspikes    = nSpikes;    % also cross-unit purposes
stat.label      = sts.label(unitsel); 
stat.freq       = sts.freq(freqindx);
stat.lfplabel   = sts.lfplabel(chansel);
stat.dimord     = 'time_lfpchan_freq';

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous sts
ft_postamble history stat


function [P] = rayleightest(x)

n = sum(~isnan(x),1);
R = resultantlength(x);    
Z = n.*R.^2;
    
P = exp(-Z).*...
(1 + (2*Z - Z.^2)./(4*n) -(24.*Z - 123*(Z.^2) + 76*(Z.^3) - 9*(Z.^4))./(288*n.^2)); %Mardia 1972

function [resLen] = resultantlength(angles)

n = sum(~isnan(angles),1);
resLen = abs(nansum(angles,1))./n;%calculate the circular variance

function [y] = ppc(crss)

dim = 1;
dof = sum(~isnan(crss),dim);
sinSum = abs(nansum(imag(crss),dim));
cosSum = nansum(real(crss),dim);
y = (cosSum.^2+sinSum.^2 - dof)./(dof.*(dof-1));     
 
function [angMean] = angularmean(angles)

angMean = angle(nansum(angles,1)); %get the mean angle


function [cfg] = trialselection(cfg,spike)

% get the number of trials or change DATA according to cfg.trials
nTrials = size(spike.trialtime,1);
if  strcmp(cfg.trials,'all')
  cfg.trials = 1:nTrials;
elseif islogical(cfg.trials)
  cfg.trials = find(cfg.trials);
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>nTrials, warning('maximum trial number in cfg.trials should not exceed length of DATA.trial')
end
if isempty(cfg.trials), error('No trials were selected');
end







