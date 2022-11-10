function [freq] = ft_spiketriggeredspectrum_stat(cfg, spike)

% FT_SPIKETRIGGEREDSPECTRUM_STAT computes phase-locking statistics for spike-LFP
% phases. These contain the PPC statistics according to Vinck et al. 2010 (Neuroimage)
% and Vinck et al. 2011 (Journal of Computational Neuroscience).
%
% Use as:
%   [stat] = ft_spiketriggeredspectrum_stat(cfg, spike)
%
% The input SPIKE should be a structure as obtained from the FT_SPIKETRIGGEREDSPECTRUM function.
%
% Configurations (cfg) 
%
% cfg.method  = string, indicating which statistic to compute. Can be:
%     -'plv' : phase-locking value, computes the resultant length over spike
%              phases. More spikes -> lower value (bias).
%     -'ang' : computes the angular mean of the spike phases.
%     -'ral' : computes the rayleigh p-value.
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
%   - double: indicates we compute statistic with a
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
% cfg.spikechannel = label of ONE unit, according to FT_CHANNELSELECTION
%
% cfg.spikesel     = 'all' (default) or numerical or logical selection of spikes.
%
% cfg.foi          = 'all' or numerical vector that specifies a subset of
%   frequencies in Hz, e.g. cfg.foi = spike.freq(1:10);                                    
%
% cfg.latency      = [beg end] in sec, or 'maxperiod',  'poststim' or
%  'prestim'.  This determines the start and end of analysis window.
%
% cfg.avgoverchan  = 'weighted', 'unweighted' (default) or 'no'.
%                  This regulates averaging of fourierspectra prior to
%                  computing the statistic.
%  - 'weighted'  : we average across channels by weighting by the LFP power.
%                  This is identical to adding the raw LFP signals in time 
%                  and then taking their FFT.
%  - 'unweighted': we average across channels after normalizing for LFP power. 
%                  This is identical to normalizing LFP signals for 
%                  their power, averaging them, and then taking their FFT.
%  - 'no'        : no weighting is performed, statistic is computed for
%                  every LFP channel.
% cfg.trials       = vector of indices (e.g., 1:2:10),
%                    logical selection of trials (e.g., [1010101010]), or
%                   'all' (default)
%
% Main outputs:
%   stat.nspikes                    =  nChancmb-by-nFreqs-nTimepoints number
%                                      of spikes used to compute stat
%   stat.dimord                     = 'chan_freq_time'
%   stat.labelcmb                   =  nChancmbs cell-array with spike vs
%                                      LFP labels
%   stat.(cfg.method)               =  nChancmb-by-nFreqs-nTimepoints  statistic
%   stat.freq                       =  1xnFreqs array of frequencies
%   stat.nspikes                    =  number of spikes used to compute
%
% The output stat structure can be plotted using ft_singleplotTFR or ft_multiplotTFR.

% Copyright (C) 2012, Martin Vinck
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


% check if the data is of spike format, and convert from old format if required
spike = ft_checkdata(spike,'datatype', 'spike', 'feedback', 'yes');

% get the options
cfg.method         = ft_getopt(cfg, 'method', 'ppc1');
cfg.channel        = ft_getopt(cfg, 'channel', 'all');
cfg.spikechannel   = ft_getopt(cfg, 'spikechannel', spike.label{1});
cfg.latency        = ft_getopt(cfg, 'latency', 'maxperiod');
cfg.spikesel       = ft_getopt(cfg, 'spikesel', 'all');
cfg.avgoverchan    = ft_getopt(cfg, 'avgoverchan', 'no');
cfg.foi            = ft_getopt(cfg, 'foi', 'all');
cfg.trials         = ft_getopt(cfg, 'trials', 'all');
cfg.timwin         = ft_getopt(cfg, 'timwin', 'all');
cfg.winstepsize    = ft_getopt(cfg, 'winstepsize', 0.01);

% ensure that the options are valid
cfg = ft_checkopt(cfg, 'method', 'char', {'ppc0', 'ppc1', 'ppc2', 'ang', 'ral', 'plv'});
cfg = ft_checkopt(cfg, 'foi',{'char', 'double'});
cfg = ft_checkopt(cfg, 'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg, 'channel', {'cell', 'char', 'double'});
cfg = ft_checkopt(cfg, 'spikesel', {'char', 'logical', 'double'});
cfg = ft_checkopt(cfg, 'avgoverchan', 'char', {'weighted', 'unweighted', 'no'});
cfg = ft_checkopt(cfg, 'latency', {'char', 'doublevector'});
cfg = ft_checkopt(cfg, 'trials', {'char', 'double', 'logical'}); 
cfg = ft_checkopt(cfg, 'timwin', {'double', 'char'}); 
cfg = ft_checkopt(cfg, 'winstepsize', {'double'}); 

cfg = ft_checkconfig(cfg, 'allowed', {'method', 'channel', 'spikechannel', 'latency', 'spikesel', 'avgoverchan', 'foi', 'trials', 'timwin', 'winstepsize'});

% collect channel information
cfg.channel        = ft_channelselection(cfg.channel, spike.lfplabel);
chansel            = match_str(spike.lfplabel, cfg.channel); 

% get the spikechannels
spikelabel       = spike.label;
cfg.spikechannel = ft_channelselection(cfg.spikechannel, spikelabel);
unitsel          = match_str(spikelabel, cfg.spikechannel);
nspikesel        = length(unitsel); % number of spike channels
if nspikesel>1, error('only one unit should be selected for now'); end
if nspikesel==0, error('no unit was selected'); end

% collect frequency information
if strcmp(cfg.foi, 'all'),  
  cfg.foi           = spike.freq(:)'; 
  freqindx          = 1:length(spike.freq);
else
  freqindx  = zeros(1,length(foi));
  for iFreq = 1:length(cfg.foi)
    freqindx(iFreq)         = nearest(spike.freq,cfg.foi(iFreq)); 
  end
end
if length(freqindx)~=length(unique(freqindx)) 
  error('Please select every frequency only once, are you sure you selected in Hz?')
end
nFreqs     = length(freqindx);
cfg.foi    = spike.freq(freqindx); % update the information again

% create the spike selection 
nSpikes = length(spike.trial{unitsel});
if strcmp(cfg.spikesel,'all'), 
  cfg.spikesel = 1:length(spike.trial{unitsel});
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
   cfg.latency = [min(spike.trialtime(:)) max(spike.trialtime(:))];
elseif strcmp(cfg.latency,'prestim')
   cfg.latency = [min(spike.trialtime(:)) 0];
elseif strcmp(cfg.latency,'poststim')
   cfg.latency = [0 max(spike.trialtime(:))];
elseif ~isrealvec(cfg.latency) && length(cfg.latency)~=2 
   error('cfg.latency should be "maxperiod", "prestim", "poststim" or 1-by-2 numerical vector');
end
if cfg.latency(1)>=cfg.latency(2), 
   error('cfg.latency should be a vector in ascending order, i.e., cfg.latency(2)>cfg.latency(1)');
end
inWindow = find(spike.time{unitsel}>=cfg.latency(1) & cfg.latency(2)>=spike.time{unitsel});

% selection of the trials
cfg        = trialselection(cfg,spike);

% do the final selection, and select on spike structure
isintrial    = ismember(spike.trial{unitsel}, cfg.trials);
spikesel     = intersect(cfg.spikesel(:),inWindow(:));
spikesel     = intersect(spikesel,find(isintrial));
cfg.spikesel = spikesel; 
spikenum     = length(cfg.spikesel); % number of spikes that were finally selected
if isempty(spikenum), warning('No spikes were selected after applying cfg.latency, cfg.spikesel and cfg.trials'); end
spike.fourierspctrm = spike.fourierspctrm{unitsel}(cfg.spikesel,chansel,freqindx);
spike.time          = spike.time{unitsel}(cfg.spikesel);
spike.trial         = spike.trial{unitsel}(cfg.spikesel);

% average the lfp channels (weighted, unweighted, or not)
if strcmp(cfg.avgoverchan,'unweighted')
  spike.fourierspctrm = spike.fourierspctrm ./ abs(spike.fourierspctrm); % normalize the angles before averaging   
  spike.fourierspctrm = nanmean(spike.fourierspctrm,2); % now rep x 1 x freq
  nChans = 1;
  outlabels = {'avgLFP'}; 
elseif strcmp(cfg.avgoverchan,'no')
  nChans = length(chansel);
  outlabels = cfg.channel;
elseif strcmp(cfg.avgoverchan,'weighted')
  spike.fourierspctrm = nanmean(spike.fourierspctrm,2); % now rep x 1 x freq
  outlabels = {'avgLFP'}; 
  nChans = 1;
end

% normalize the spectrum first
spike.fourierspctrm = spike.fourierspctrm ./ abs(spike.fourierspctrm); % normalize the angles before averaging   
ft_progress('init', 'text',     'Please wait...');
if strcmp(cfg.timwin,'all')
  freq.time = 'all';
  switch cfg.method
    case 'ang'
     out  = angularmean(spike.fourierspctrm);
    case 'plv'
     out  = resultantlength(spike.fourierspctrm);
    case 'ral'
     out  = rayleightest(spike.fourierspctrm); % the rayleigh test
    case 'ppc0'
     out  = ppc(spike.fourierspctrm);
    case {'ppc1', 'ppc2'}
      
      % check the final set of trials present in the spikes
      trials = unique(spike.trial);

      % loop init for PPC 2.0
      [S,SS,dof,dofSS] = deal(zeros(1,nChans,nFreqs));
      nTrials = length(trials);
      if nTrials==1
        warning('computing ppc1 or ppc2 can only be performed with more than 1 trial');
      end
      for iTrial = 1:nTrials % compute the firing rate
          trialNum      = trials(iTrial);
          spikesInTrial = find(spike.trial == trialNum);
          spc           = spike.fourierspctrm(spikesInTrial,:,:);
          ft_progress(iTrial/nTrials, 'Processing trial %d from %d', iTrial, nTrials);              
          % compute PPC 1.0 and 2.0 according to Vinck et al. (2011) using summation per trial
          if strcmp(cfg.method,'ppc2')
            if ~isempty(spc)
              m = nanmean(spc,1); % no problem with NaN
              hasNum = ~isnan(m);
              S(hasNum)  = S(hasNum)  + m(hasNum); % add the imaginary and real components
              SS(hasNum) = SS(hasNum) + m(hasNum).*conj(m(hasNum));
              dof(hasNum) = dof(hasNum) + 1; % dof needs to be kept per frequency
            end
          elseif strcmp(cfg.method,'ppc1')
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
      elseif strcmp(cfg.method, 'ppc2')
        out(hasTrl) = (S(hasTrl).*conj(S(hasTrl)) - SS(hasTrl))./(dof(hasTrl).*(dof(hasTrl)-1));
      end
  end
  nSpikes = sum(~isnan(spike.fourierspctrm));
else % compute time-resolved spectra of statistic
  
  % make the sampling axis for the window
  bins      = cfg.latency(1):cfg.winstepsize:cfg.latency(2);
  N         = length(bins)-1; % number of bins
  wintime   = 0:cfg.winstepsize:cfg.timwin;
  win       = ones(1,length(wintime));
  freq.time = (bins(2:end)+bins(1:end-1))/2;
  if ~mod(length(win),2), win = [win 1]; end   % make sure the number of samples is uneven.

  out      = NaN(N,nChans,nFreqs);
  nSpikes  = zeros(N,nChans,nFreqs);
  for iChan = 1:nChans
    for iFreq = 1:nFreqs
  
      spctra       = spike.fourierspctrm(:,iChan,iFreq); % values to accumulate
      tm           = spike.time;
      hasnan       = isnan(spctra);    
      tm(hasnan) = [];
      spctra(hasnan) = [];

      % compute the degree of freedom per time bin and the index for every spike
      [dof, indx] = histc(tm, bins); % get the index per spike, and number per bin                          
      if isempty(dof), continue,end    
      toDel       = indx==length(bins) | indx==0; % delete those points that are equal to last output histc or don't fall in
      spctra(toDel) = [];
      indx(toDel)   = [];
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
          trials  = unique(spike.trial);
          nTrials = length(trials);
          [S,SS,dofS,dofSS] = deal(zeros(length(bins)-1,1));
          if nTrials==1
             warning('computing ppc1 or ppc2 can only be performed with more than 1 trial');
          end      
          % compute the new ppc versions
          for iTrial = 1:nTrials

            ft_progress(iTrial/nTrials, 'Processing trial %d from %d for freq %d and chan %d', iTrial, nTrials, iFreq, iChan);                
             % select the spectra, time points, and trial numbers again
             trialNum      = trials(iTrial);
             spikesInTrial = find(spike.trial == trialNum);
             if isempty(spikesInTrial), continue,end
             spctraTrial  = spike.fourierspctrm(spikesInTrial,iChan,iFreq);       
             tm           = spike.time;
             hasnan       = isnan(spctraTrial);    
             tm(hasnan) = [];
             spctraTrial(hasnan) = [];

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
ft_progress('close');

% collect the outputs: in labelcmb representation
outparam        = cfg.method;
freq.(outparam) = permute(out,[2 3 1]);
freq.nspikes    = permute(nSpikes,[2 3 1]);    % also cross-unit purposes
freq.labelcmb = cell(nChans,2);
freq.labelcmb(1:nChans,1) = cfg.spikechannel;
for iCmb = 1:nChans
  freq.labelcmb{iCmb,2}   = outlabels{iCmb}; 
end  
freq.freq       = spike.freq(freqindx);
freq.dimord     = 'chancmb_freq_time';

% do the general cleanup and bookkeeping at the end of the function

ft_postamble previous   spike
ft_postamble provenance freq
ft_postamble history    freq


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P] = rayleightest(x)

n = sum(~isnan(x),1);
R = resultantlength(x);    
Z = n.*R.^2;
    
P = exp(-Z).*...
(1 + (2*Z - Z.^2)./(4*n) -(24.*Z - 123*(Z.^2) + 76*(Z.^3) - 9*(Z.^4))./(288*n.^2)); %Mardia 1972

function [resLen] = resultantlength(angles)

n = sum(~isnan(angles),1);
resLen = abs(nansum(angles,1))./n; %calculate the circular variance

function [y] = ppc(crss)

dim = 1;
dof = sum(~isnan(crss),dim);
sinSum = abs(nansum(imag(crss),dim));
cosSum = nansum(real(crss),dim);
y = (cosSum.^2+sinSum.^2 - dof)./(dof.*(dof-1));     
 
function [angMean] = angularmean(angles)

angMean = angle(nanmean(angles,1)); %get the mean angle


function [cfg] = trialselection(cfg,spike)

% get the number of trials or change DATA according to cfg.trials
nTrials = size(spike.trialtime,1);
if  strcmp(cfg.trials, 'all')
  cfg.trials = 1:nTrials;
elseif islogical(cfg.trials) || all(cfg.trials==0 | cfg.trials==1)
  cfg.trials = find(cfg.trials);
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>nTrials, warning('maximum trial number in cfg.trials should not exceed length of DATA.trial')
end
if isempty(cfg.trials), error('No trials were selected');
end

function m = nansum(x,dim)

% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    n = sum(~nans);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x);
else
    % Count up non-NaNs.
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x,dim);
end






