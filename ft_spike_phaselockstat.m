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
%   STS should be a structure as obtained from from the FT_SPIKE_TRIGGEREDSPECTRUM function. 
%
%   Configurations (CFG) related to calculation of phaselockingvalue
%
%   cfg.channel                     = Nx1 cell-array or numerical array with selection of 
%                                     channels (default = 'all'),See CHANNELSELECTION for details
%   cfg.spikesel                    = 'all' (default) or numerical or logical selection of spikes.
%   cfg.foi                         = 'all' or numerical vector that specifies a subset of
%                                     frequencies in Hz (and not in indices!).
%   cfg.latency                     = [beg end] in sec, or 'maxperiod', 'poststim' or 'prestim'.
%                                     This determines the start and end of time-vector.
%   cfg.spikechannel                = label of ONE unit, according to FT_CHANNELSELECTION
%   cfg.chanavg                     = 'yes' or 'no' (default).
%   cfg.powweighted                 = 'yes' or 'no' (default). If 'yes', we
%                                     average across channels by weighting by the LFP power.
%   Main outputs:
%
%     stat.ppc0                       =  nChan-by-nFreqs matrix with the ppc 0
%     stat.ppc1                       =  nChan-by-nFreqs matrix with the ppc 1
%     stat.ppc2                       =  nChan-by-nFreqs matrix with the ppc 2
%     stat.ral                        =  nChan-by-nFreqs rayleigh stat.
%     stat.ang                        =  nChan-by-nFreqs mean phase angle.
%     stat.dofspike                   =  nChan-by-nFreqs number of spikes
%     stat.spikechannel               =  name of unit;
%     stat.plv                        =  nChan-by-nFreqs phase-locking value
%     
%   References:
%   Vinck et al. (2010) Neuroimage
%   Vinck et al. (2011) Journal of Computational Neuroscience

%   Copyright (c) Martin Vinck (2011), University of Amsterdam.
%
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'foi','t_ftimwin'});
if  isequal(cfg.taper, 'dpss')
  cfg = ft_checkconfig(cfg, 'required', {'tapsmofrq'});
  cfg = ft_checkopt(cfg,'tapsmofrq',{'doublevector', 'doublescalar'});
end

% get the options
cfg.channel        = ft_getopt(cfg,'channel', 'all');
cfg.spikechannel   = ft_getopt(cfg,'spikechannel', sts.spikechannel{1});
cfg.latency        = ft_getopt(cfg,'latency', 'maxperiod');
cfg.spikesel       = ft_getopt(cfg,'spikesel', 'all');
cfg.chanavg        = ft_getopt(cfg,'chanavg', 'no');
cfg.powweighted    = ft_getopt(cfg,'powweighted', 'no');
cfg.foi            = ft_getopt(cfg,'foi', 'all');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'foi',{'char', 'double'});
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'channel', {'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'spikesel', {'char', 'logical', 'double'});
cfg = ft_checkopt(cfg,'chanavg', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'powweighted', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'latency', {'char', 'doublevector'});

% collect channel information
cfg.channel        = ft_channelselection(cfg.channel, sts.label);
chansel            = match_str(sts.label, cfg.channel); 

% get the spikechannels
spikelabel       = sts.spikechannel;
cfg.spikechannel = ft_channelselection(cfg.spikechannel, spikelabel);
unitsel          = match_str(spikelabel, cfg.spikechannel);
nspikesel        = length(unitsel); % number of spike channels
if nspikesel>1, error('only one unit should be selected for now'); end

% collect frequency information
if strcmp(cfg.foi, 'all'),  
  cfg.foi           = sts.freq; 
elseif ~isrealvec(cfg.foi), 
  error('MATLAB:fieldtrip:spike_phaselocking:cfg:foi:unknownOption',...
  'cfg.foi should be "all" or vector of frequencies in Hz')
end
freqindx         = nearest_nd(sts.freq,cfg.foi); 
if length(freqindx)~=length(unique(freqindx)) 
  error('MATLAB:fieldtrip:spike_phaselocking:cfg:foi:notUniqueSelection',... 
  'Please select every frequency only once, are you sure you selected in Hz?')
end
nFreqs     = length(freqindx);
cfg.foi    = sts.freq(freqindx); % update the information again

% create the spike selection (trial selection can be done via this option)
nSpikes = length(sts.trial{unitsel});
if strcmp(cfg.spikesel,'all'), 
  cfg.spikesel = 1:length(sts.trial{unitsel});
elseif islogical(cfg.spikesel)
  cfg.spikesel = find(cfg.spikesel);
elseif ~isrealvec(cfg.spikesel) 
  error('MATLAB:fieldtrip:spike_phaselocking:cfg:spikesel',... 
  'cfg.spikesel should be a numerical or logical vector ')
end
if ~isempty(cfg.spikesel)
if max(cfg.spikesel)>nSpikes || length(cfg.spikesel)>nSpikes % ease debugging of this error
  error('MATLAB:fieldtrip:spike_phaselocking:cfg:spikesel',...
  'cfg.spikesel must not exceed number of spikes and select every spike just once')
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
   error('MATLAB:fieldtrip:spike_phaselocking:cfg:latency',...
    'cfg.latency should be "maxperiod", "prestim", "poststim" or 1-by-2 numerical vector');
end
if cfg.latency(1)>=cfg.latency(2), 
   error('MATLAB:fieldtrip:spike_phaselocking:cfg:latency:wrongOrder',...
   'cfg.latency should be a vector in ascending order, i.e., cfg.latency(2)>cfg.latency(1)');
end
inWindow = find(sts.time{unitsel}>=cfg.latency(1) & cfg.latency(2)>=sts.time{unitsel});

% do the final selection
%isNum        = find(isfinite(sts.fourierspctrm{unitsel}(:,1,1)));
spikesel     = intersect(cfg.spikesel(:),inWindow(:));
cfg.spikesel = spikesel; %intersect(spikesel(:),isNum(:));
spikenum     = length(cfg.spikesel); % number of spikes that were finally selected
if isempty(spikenum), warning('MATLAB:fieldtrip:spike_phaselocking:silentNeuron',...
  'No spikes were selected after applying cfg.latency, cfg.spikesel and cfg.trials');
end
sts.fourierspctrm{unitsel} = sts.fourierspctrm{unitsel}(cfg.spikesel,chansel,freqindx);
sts.time{unitsel}          = sts.time{unitsel}(cfg.spikesel);
sts.trial{unitsel}         = sts.trial{unitsel}(cfg.spikesel);

if strcmp(cfg.powweighted,'no')
  sts.fourierspctrm{unitsel} = sts.fourierspctrm{unitsel} ./ abs(sts.fourierspctrm{unitsel}); % normalize the angles
  if strcmp(cfg.chanavg,'yes')
    sts.fourierspctrm{unitsel} = nanmean(sts.fourierspctrm{unitsel},2); % now rep x 1 x freq
%    sts.fourierspctrm{unitsel} = sts.fourierspctrm{unitsel} ./ abs(sts.fourierspctrm{unitsel}); % normalize the angles  
    nChans = 1;
  else
    nChans = length(chansel);
  end
else
  if strcmp(cfg.chanavg,'yes')
    sts.fourierspctrm{unitsel} = nanmean(sts.fourierspctrm{unitsel},2); % now rep x 1 x freq
    nChans = 1;
  else
    nChans = length(chansel);
  end
end 
  
% implement the new function
trials = unique(sts.trial{unitsel});

% loop init for PPC 2.0
[S,SS,dof,S1,SS1,dofS1,dofSS1,S2,SS2,dofS2,dofSS2] = deal(zeros(1,nChans,nFreqs)); 

nTrials = length(trials);
for iTrial = 1:nTrials % compute the firing rate
    trialNum      = trials(iTrial);
    spikesInTrial = find(sts.trial{unitsel} == trialNum);
    spcU           = sts.fourierspctrm{unitsel}(spikesInTrial,:,:);
    spc           = spcU./abs(spcU);
    
    % compute PPC 2.0
    if ~isempty(spc)
      m = nanmean(spc,1); % no problem with NaN
      hasNum = ~isnan(m);
      S(hasNum)  = S(hasNum)  + m(hasNum); % add the imaginary and real components
      SS(hasNum) = SS(hasNum) + m(hasNum).*conj(m(hasNum));
      dof(hasNum) = dof(hasNum) + 1; % dof needs to be kept per frequency
    end
    
    % compute PPC 1.0
    if ~isempty(spc)
       n      = sum(~isnan(spc),1);
       m      = nansum(spc,1); 
       hasNum = ~isnan(m);    
       S1(hasNum)  = S1(hasNum)  + m(hasNum); % add the imaginary and real components
       SS1(hasNum) = SS1(hasNum) + m(hasNum).*conj(m(hasNum));
       dofS1(hasNum)  = dofS1(hasNum)  + n(hasNum);
       dofSS1(hasNum) = dofSS1(hasNum) + n(hasNum).^2;
    end                              
end

[ppc1,ppc2,ppc3]   = deal(NaN(1,nChans,nFreqs));
hasTrl = dof>1;
ppc1(hasTrl) = (S1(hasTrl).*conj(S1(hasTrl)) - SS1(hasTrl))./(dofS1(hasTrl).^2 - dofSS1(hasTrl));
ppc2(hasTrl) = (S(hasTrl).*conj(S(hasTrl)) - SS(hasTrl))./(dof(hasTrl).*(dof(hasTrl)-1));

ral  = rayleightest(sts.fourierspctrm{unitsel},1); % the rayleigh test
ppc0 = ppc(angle(sts.fourierspctrm{unitsel}),1); % the old ppc
ang  = angularmean(sts.fourierspctrm{unitsel},1);

% already compute the quantities that are used for cross-unit PPC 
stat.doftrial = dof;
dof    = sum(~isnan(sts.fourierspctrm{unitsel}),1);

stat.ppc0   = ppc0;
stat.ppc1   = ppc1;
stat.ppc2   = ppc2;
stat.spikechannel = spikelabel{unitsel};
stat.dofspike = dof;    % also cross-unit purposes
stat.ang    = ang;
stat.ral    = ral;
stat.plv    = resultantlength(sts.fourierspctrm{unitsel},1);
stat.freq   = sts.freq(freqindx);
stat.label  = sts.label(chansel);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous sts
ft_postamble history stat


function [P] = rayleightest(x,dim)

%RAYLEIGHTEST   Rayleigh test
%   
%   P = RAYLEIGHTEST(X,DIM) takes the rayleightest along the dimension DIM of X. 
%
%   Example: If X = 2*pi*rand(100,5,10);
%
%   then p = rayleightest(X,1) returns a matrix of probabilities of size
%   (5,10) under the null hypothesis that the data is drawn from a uniform
%   distribution
%
%   Class support for input X:
%      float: double (both complex and real)
%
%   See also RESULTANTLENGTH
%   Copyright 2008 Martin Vinck

n = sum(~isnan(x),dim);
R = resultantlength(x,dim);    
Z = n.*R.^2;
    
P = exp(-Z).*...
(1 + (2*Z - Z.^2)./(4*n) -(24.*Z - 123*(Z.^2) + 76*(Z.^3) - 9*(Z.^4))./(288*n.^2)); %Mardia 1972

function [resLen] = resultantlength(angles,dim)

if nargin==1, 
  dim = min(find(size(angles)~=1));
  if isempty(dim), dim = 1; end
end  

if isreal(angles)
    angles = exp(1i*angles);
else
    angles = angles./abs(angles);
end

n = sum(~isnan(angles),dim);

resLen = abs(nansum(angles,dim))./n;%calculate the circular variance


function [y] = ppc(angles,dim)

if nargin<2
  dim = 1;
end

crss = exp(1i*angles);
dof = sum(~isnan(crss),dim);
sinSum = abs(nansum(imag(crss),dim));
cosSum = nansum(real(crss),dim);
y = (cosSum.^2+sinSum.^2 - dof)./(dof.*(dof-1));     
 

function [angMean] = angularmean(angles,dim)

% works with NaN as well...
if nargin==1, 
  dim = min(find(size(angles)~=1));
  if isempty(dim), dim = 1; end
end  
    
if isreal(angles)
    angles = exp(i*angles);
else
    angles = angles./abs(angles);
end
angMean = angle(nansum(angles,dim)); %get the mean angle






