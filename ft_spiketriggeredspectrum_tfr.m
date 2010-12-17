function [sts_tfr] = ft_spiketriggeredspectrum_ppc_tfr(cfg,freq)

%   FT_SPIKETRIGGEREDSPECTRUM_TFR computes time-frequency representation of PPC, phase and rayleigh test
%   Getting a TFR from spike phases is complicated, because spike numbers may strongly vary over
%   time. However, this is solved with the PPC statistic (Vinck, 2010, Neuroimage), being unbiased
%   by the number of spikes.
% 
% Inputs:
%   STS should be a structure as obtained from from the FT_SPIKETRIGGEREDSPECTRUM function. 
%
%   Configurations (CFG)
%
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
%                                     
%   Main outputs:
%
%     sts_tfr.ppc                =  nChan-by-nFreqs matrix with the raw ppc spectrum.
%     sts_tfr.phase              =  nChan-by-nFreqs matrix with the raw phases
%     sts_tfr.ral                =  nChan-by-nFreqs matrix with the rayleigh test p-values
%     sts_tfr.time               =  vector with center timepoints, i.e., the windows are
%                                   centered on these timepoints
%   Other outputs as usual.
%
%   See also FT_SPIKETRIGGEREDSPECTRUM
%

%   Copyright (c) Martin Vinck
% NOTE: I left cfg.trials away, since cfg.spikesel can do what cfg.trials can do but in a more
% general way. getting cfg.spikesel from a given set of trials is rather easy - so I think we should
% let users do this themselves

% defaults business
defaults.channel      = {'all'};             
defaults.latency      = {'maxperiod'};              
defaults.spikesel     = {'all'};
defaults.foilim       = {'all'};
defaults.winlen       = {0.1};
defaults.fsample      = {1000};
cfg = ft_spike_sub_defaultcfg(cfg,defaults);

% channel selection business
cfg.channel        = ft_channelselection(cfg.channel, freq.label);
chansel            = match_str(freq.label, cfg.channel); 

% frequency selection business
if strcmp(cfg.foilim, 'all'),  
  cfg.foilim           = [freq.freq(1) freq.freq(end)]; 
elseif ~isrealvec(cfg.foilim) 
  error('MATLAB:fieldtrip:spike_phaselocking:cfg:foi:unknownOption',...
  'cfg.foi should be "all" or vector of frequencies in Hz')
end
fbeg             = nearest(freq.freq,cfg.foilim(1));
fend             = nearest(freq.freq,cfg.foilim(end));
freqindx         = fbeg:fend;
freqindx         = unique(freqindx);
cfg.foilim       = freq.freq(freqindx); % update the information again
nFreqs           = length(cfg.foilim);

% create the spike selection (trial selection can be done via this option
% we need to re-enter the trial selection again!
nSpikes = length(freq.trial);
if strcmp(cfg.spikesel,'all'), 
  cfg.spikesel = 1:length(freq.trial);
elseif islogical(cfg.spikesel)
  cfg.spikesel = find(cfg.spikesel);
elseif ~isrealvec(cfg.spikesel) 
  error('MATLAB:fieldtrip:spike_phaselocking:cfg:spikesel',... 
  'cfg.spikesel should be a numerical or logical vector ')
end
if max(cfg.spikesel)>nSpikes || length(cfg.spikesel)>nSpikes % ease debugging of this error
  error('MATLAB:fieldtrip:spike_phaselocking:cfg:spikesel',...
  'cfg.spikesel must not exceed number of spikes and select every spike just once')
end
  
% select on basis of latency, we use the new format, with freq.trialtime in it
begTrialLatency = freq.trialtime; 
endTrialLatency = freq.trialtime;  
if strcmp(cfg.latency, 'maxperiod'),
   cfg.latency = [min(begTrialLatency) max(endTrialLatency)];
elseif strcmp(cfg.latency,'minperiod')
   cfg.latency = [max(begTrialLatency) min(endTrialLatency)];
elseif strcmp(cfg.latency,'prestim')
   cfg.latency = [min(begTrialLatency) 0]; % ERROR HERE, CHECK OTHER FUNCS!
elseif strcmp(cfg.latency,'poststim')
   cfg.latency = [0 max(endTrialLatency)];
elseif ~isrealvec(cfg.latency) && length(cfg.latency)~=2 
   error('MATLAB:fieldtrip:spike_phaselocking:cfg:latency',...
    'cfg.latency should be "maxperiod", "prestim", "poststim" or 1-by-2 numerical vector');
end
if cfg.latency(1)>=cfg.latency(2), 
   error('MATLAB:fieldtrip:spike_phaselocking:cfg:latency:wrongOrder',...
   'cfg.latency should be a vector in ascending order, i.e., cfg.latency(2)>cfg.latency(1)');
end
inWindow = find(freq.time>=cfg.latency(1) & cfg.latency(2)>=freq.time);

% do the final selection
cfg.spikesel = intersect(cfg.spikesel(:),inWindow(:));
spikenum     = length(cfg.spikesel); % number of spikes that were finally selected
if isempty(spikenum), warning('MATLAB:fieldtrip:spike_phaselocking:silentNeuron',...
  'No spikes were selected after applying cfg.latency, cfg.spikesel and cfg.trials');
end
freq.fourierspctrm = freq.fourierspctrm(cfg.spikesel,chansel,freqindx);
freq.time          = freq.time(cfg.spikesel);
freq.trial         = freq.trial(cfg.spikesel);

% normalize the fourier spectrum to be unity complex vectors
freq.fourierspctrm = freq.fourierspctrm ./ abs(freq.fourierspctrm); 

% take care of all the time business and window business
dt      = 1/(cfg.fsample);
bins    = cfg.latency(1):dt:cfg.latency(2);
N       = length(bins)-1; % number of bins
wintime = 0:dt:cfg.winlen;
win     = ones(1,length(wintime));

% make sure the number of samples is uneven.
if ~mod(length(win),2)
  win = [win 1];
end

% find the indices belonging to the spikes
[dof, indx] = histc(freq.time, bins); % get the index per spike, and number per bin                          
dof(end) = []; % the last bin is a single point in time, so we delete it
dof = dof(:); % force it to be row
dof  = conv(dof,win,'same'); % get the total number of spikes across all trials
toDel = indx==length(bins) | indx==0; % delete those points that are equal to last output histc
freq.fourierspctrm(toDel,:,:) = []; % delete those spikes from fourierspctrm as well
indx(toDel) = []; % make sure index doesn't contain them

% preallocate the ppc and ral matrix
nChans  = length(chansel);
nFreqs  = length(freqindx);
[ppc,phs,ral] = deal(zeros(N,nChans,nFreqs));
for iChan = 1:nChans
  for iFreq = 1:nFreqs
    
    % we will simply sum all the complex numbers across time / spikes
    vals = freq.fourierspctrm(:,chansel(iChan),freqindx(iFreq)); % values to accumulate
    vals(~isnan(vals)) = 0; % set NaNs to zero, that way they won't affect the sum    
    x    = accumarray(indx(:),vals,[N 1]); % simply the sum of the complex numbers
    y    = conv(x,win,'same'); % convolution again just means a sum                                
    
    % compute the ppc using a new trick, doesn't require any loops
    ppc(:,iChan,iFreq) = (y.*conj(y) - dof)./(dof.*(dof-1)); % simplest form of ppc
    phs(:,iChan,iFreq) = angle(y); % 
    Z = abs(y).^2./dof;
    P = exp(-Z).*...
    (1 + (2*Z - Z.^2)./(4*dof) -(24.*Z - 123*(Z.^2) + 76*(Z.^3) - 9*(Z.^4))./(288*dof.^2)); %Mardia 1972
    ral(:,iChan,iFreq) = P;     
  end
end

% gather the results
sts_tfr.ppc   = ppc;
sts_tfr.phase = phs;
sts_tfr.ral   = ral;
sts_tfr.freq  = freq.freq(freqindx);
sts_tfr.label = freq.label(chansel);
sts_tfr.time  = bins(1:end-1) + 0.5*dt; % center time-points
sts_tfr.cfg   = cfg;
    
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
if isfield(freq,'cfg'),cfg.previous = freq.cfg; end
% remember the exact configuration details in the output 
sts_tfr .cfg     = cfg;
