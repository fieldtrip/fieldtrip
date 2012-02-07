function [sts_tfr] = ft_spiketriggeredspectrum_tfr(cfg, sts)

% FT_SPIKETRIGGEREDSPECTRUM_TFR computes time-frequency representation of PPC,
% phase and rayleigh test Getting a TFR from spike phases is complicated,
% because spike numbers may strongly vary over time. However, this is solved
% with the PPC statistic (Vinck et al, 2010; Neuroimage), Vinck et al., 2011,
% Journal of Computational Neuroscience, being unbiased by the number of
% spikes.
%
% Use as
%   [stat] = ft_spiketriggeredspectrum_ppc_tfr(cfg,stat)
%
% Inputs:
%   STS should be a structure as obtained from from the
%   FT_SPIKETRIGGEREDSPECTRUM or FT_SPIKE_TRIGGEREDSPECTRUM function
%
% Configurations:
%   cfg.channel                     = Nx1 cell-array or numerical array with selection of
%                                     channels (default = 'all'). See FT_CHANNELSELECTION for details
%   cfg.spikesel                    = 'all' (default) or numerical or logical selection of spikes.
%   cfg.foilim                      = [fbeg fend] in Hz.
%   cfg.latency                     = [beg end] in sec, or 'maxperiod' (default), 'poststim' or 'prestim'.
%                                     This determines the start and end of time-vector.
%   cfg.fsample                     = sampling frequency of time-axis:
%                                     e.g., cfg.fsample = 1000 will make time-points to be separated
%                                     with 1/1000 = 0.001 sec.
%   cfg.winlen                      = length of the window in which we compute ppc/phase in seconds.
%   cfg.chanavg                     = 'yes' or 'no' (default): average over chans.
%   cfg.spikechannel                = string or cell of single channel to
%                                     compute stats for
%
% See also FT_SPIKETRIGGEREDSPECTRUM, FT_SPIKE_PHASELOCKSTAT,
% FT_SPIKE_TRIGGEREDSPECTRUM

% Copyright (C) 2010, Martin Vinck
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
cfg.channel        = ft_getopt(cfg,'channel', 'all');
cfg.spikechannel   = ft_getopt(cfg,'spikechannel', sts.label{1});
cfg.latency        = ft_getopt(cfg,'latency', 'maxperiod');
cfg.spikesel       = ft_getopt(cfg,'spikesel', 'all');
cfg.chanavg        = ft_getopt(cfg,'chanavg', 'no');
cfg.foilim         = ft_getopt(cfg,'foilim', 'all');
cfg.fsample        = ft_getopt(cfg,'fsample', 1000);
cfg.winlen         = ft_getopt(cfg,'winlen', 0.1);

cfg = ft_checkopt(cfg,'channel', {'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double', 'empty'});
cfg = ft_checkopt(cfg,'winlen', 'double');
cfg = ft_checkopt(cfg,'latency', {'double', 'char'});
cfg = ft_checkopt(cfg,'chanavg', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'fsample', 'double');
cfg = ft_checkopt(cfg,'foilim', {'char', 'doublevector'});
cfg = ft_checkopt(cfg,'spikesel', {'logical', 'double', 'char'});

% get the spikechannels
spikelabel       = sts.label;
cfg.spikechannel = ft_channelselection(cfg.spikechannel, spikelabel);
unitsel          = match_str(spikelabel, cfg.spikechannel);
nspikesel        = length(cfg.spikechannel); % number of spike channels
if nspikesel>1, error('only one unit should be selected for now'); end

sts.trial = sts.trial{unitsel};
sts.time  = sts.time{unitsel};
sts.fourierspctrm = sts.fourierspctrm{unitsel};

% channel selection business
cfg.channel        = ft_channelselection(cfg.channel, sts.lfplabel);
chansel            = match_str(sts.lfplabel, cfg.channel); 

% frequency selection business
if strcmp(cfg.foilim, 'all'),  
  cfg.foilim           = [sts.freq(1) sts.freq(end)]; 
end
fbeg             = nearest(sts.freq,cfg.foilim(1));
fend             = nearest(sts.freq,cfg.foilim(end));
freqindx         = fbeg:fend;
freqindx         = unique(freqindx);
cfg.foilim       = sts.freq(freqindx); % update the information again
nFreqs           = length(cfg.foilim);

% create the spike selection (trial selection can be done via this option
% we need to re-enter the trial selection again!
nSpikes = length(sts.trial);
if strcmp(cfg.spikesel,'all'), 
  cfg.spikesel = 1:length(sts.trial);
elseif islogical(cfg.spikesel)
  cfg.spikesel = find(cfg.spikesel);
end
if max(cfg.spikesel)>nSpikes | length(cfg.spikesel)>nSpikes % ease debugging of this error
  error('MATLAB:fieldtrip:spike_phaselocking:cfg:spikesel',...
  'cfg.spikesel must not exceed number of spikes and select every spike just once')
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
inWindow = find(sts.time>=cfg.latency(1) & cfg.latency(2)>=sts.time);

% do the final selection
cfg.spikesel     = intersect(cfg.spikesel(:),inWindow(:));
spikenum     = length(cfg.spikesel); % number of spikes that were finally selected
if isempty(spikenum), warning('MATLAB:fieldtrip:spike_phaselocking:silentNeuron',...
  'No spikes were selected after applying cfg.latency, cfg.spikesel and cfg.trials');
end
sts.fourierspctrm = sts.fourierspctrm(cfg.spikesel,chansel,freqindx);
sts.time      = sts.time(cfg.spikesel);
sts.trial     = sts.trial(cfg.spikesel);

% normalize the fourier spectrum to be unity complex vectors
sts.fourierspctrm = sts.fourierspctrm ./ abs(sts.fourierspctrm); 
if strcmp(cfg.chanavg,'yes')
  sts.fourierspctrm = nanmean(sts.fourierspctrm,2); % now rep x 1 x freq
  sts.fourierspctrm = sts.fourierspctrm ./ abs(sts.fourierspctrm); % normalize the angles  
  chansel            = 1;
  sts.label         = {'avgchan'};
end

% take care of all the time business and window business
dt      = 1/(cfg.fsample);
bins    = cfg.latency(1):dt:cfg.latency(2);
N    = length(bins)-1; % number of bins
wintime = 0:dt:cfg.winlen;
win     = ones(1,length(wintime));

% make sure the number of samples is uneven.
if ~mod(length(win),2)
  win = [win 1];
end

% avoid over-vectorizing the code here. MATLAB implementation should be fast enough here.
% I prefer low memory load (< 2GB ram) over a few factors of speed.
nChans  = length(chansel);
nFreqs  = length(freqindx);
[ppc0,ppc1,ppc2,phs,ral,plv] = deal(NaN(N,nChans,nFreqs));
df = zeros(N,nChans,nFreqs);
trials = unique(sts.trial);
nTrials = length(trials);
for iChan = 1:nChans
  for iFreq = 1:nFreqs
    
    % we will simply sum all the complex numbers across time / spikes
    vals   = sts.fourierspctrm(:,iChan,iFreq); % values to accumulate
    hasnan = isnan(vals);    
    vals(hasnan) = []; % since accumarray is essentially a summing procedure.
        
    tm  = sts.time;
    trl = sts.trial;
    tm(hasnan)  = [];
    trl(hasnan) = [];
    
    % find the indices belonging to the spikes
    [dof, indx] = histc(tm, bins); % get the index per spike, and number per bin                          
    if isempty(dof), continue,end    
    dof(end)    = []; % the last bin is a single point in time, so we delete it
    dof         = dof(:); % force it to be row
    dof         = conv2(dof(:),win(:),'same'); % get the total number of spikes across all trials
    toDel       = indx==length(bins) | indx==0; % delete those points that are equal to last output histc or don't fall in
    vals(toDel) = [];
    tm(toDel)   = [];
    trl(toDel)  = [];
    indx(toDel) = [];
    df(:,iChan,iFreq) = dof;
    
    x    = accumarray(indx(:),vals,[N 1]); % simply the sum of the complex numbers
    y    = conv2(x(:),win(:),'same'); % convolution again just means a sum                                
    
    % compute the ppc using a new trick, doesn't require any loops
    hasnum = dof>1;
    ppc0(hasnum,iChan,iFreq) = (y(hasnum).*conj(y(hasnum)) - dof(hasnum))./(dof(hasnum).*(dof(hasnum)-1)); % simplest form of ppc
    plv(hasnum,iChan,iFreq) = abs(y(hasnum)./dof(hasnum));
    hasnum0 = dof>0;        
    phs(hasnum0,iChan,iFreq) = angle(y(hasnum0)); % 
    Z = abs(y).^2./dof;
    P = exp(-Z).*...
    (1 + (2*Z - Z.^2)./(4*dof) -(24.*Z - 123*(Z.^2) + 76*(Z.^3) - 9*(Z.^4))./(288*dof.^2)); %Mardia 1972
    ral(hasnum,iChan,iFreq) = P(hasnum);     
    
    [S_norm,SS_norm,dofS_norm,S,SS,dofSS,dofS] = deal(zeros(length(bins)-1,1));

    % compute the new ppc versions
    for iTrial = 1:nTrials
       trialNum      = trials(iTrial);
       spikesInTrial = find(trl == trialNum);
       if ~isempty(spikesInTrial)
         trlvals = vals(spikesInTrial);       
         [d, id] = histc(tm(spikesInTrial), bins); % get the index per spike, and number per bin                          
         d(end) = [];
         toDel = id==length(bins) | id==0; % delete those points that are equal to last output histc
         trlvals(toDel,:,:) = []; % delete those spikes from fourierspctrm as well
         id(toDel) = []; % make sure index doesn't contain them
         x    = accumarray(id(:),trlvals,[N 1]); % simply the sum of the complex numbers
         y    = conv2(x(:),win(:),'same'); % convolution again just means a sum                                  
         d    = conv2(d(:),win(:),'same'); % get the dof per trial

         S     = S  + y; % add the imaginary and real components
         SS    = SS + y.*conj(y);
         dofS  = dofS  + d(:);
         dofSS = dofSS + d(:).^2;
       
         sl                  = d(:)>0;
         m                   = y./d(:);
         S_norm(sl)          = S_norm(sl)  + m(sl); % add the imaginary and real components
         SS_norm(sl)         = SS_norm(sl) + m(sl).*conj(m(sl));
         dofS_norm(sl)       = dofS_norm(sl) + 1;                     
       end
    end
    % we need at least two trials
    hasNum = dofS_norm>1;
    ppc1(hasNum,iChan,iFreq) = (S(hasNum).*conj(S(hasNum)) - SS(hasNum))./(dofS(hasNum).^2 - dofSS(hasNum));
    ppc2(hasNum,iChan,iFreq) = (S_norm(hasNum).*conj(S_norm(hasNum)) - SS_norm(hasNum))./(dofS_norm(hasNum).*(dofS_norm(hasNum)-1));
  end
end

% gather the results
sts_tfr.ppc0     = ppc0;
sts_tfr.ppc1     = ppc1;
sts_tfr.ppc2     = ppc2;
sts_tfr.plv      = plv;
sts_tfr.label    = spikelabel(unitsel);
sts_tfr.ang      = phs;
sts_tfr.ral      = ral;
sts_tfr.dof      = df;
sts_tfr.freq     = sts.freq(freqindx);
sts_tfr.lfplabel = sts.lfplabel(chansel);
sts_tfr.time     = bins(1:end-1) + 0.5*dt; % center time-points
sts_tfr.cfg      = cfg;
sts_tfr.dimord   = 'time_lfpchan_freq';

