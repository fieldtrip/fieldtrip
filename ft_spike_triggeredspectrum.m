function [freq] = ft_spike_triggeredspectrum(cfg, data, spike)

% FT_SPIKE_TRIGGEREDSPECTRUM computes the Fourier spectrup of the LFP around the spikes.
%
% When spikes have binary representation, as obtained from APPENDSPIKE or FT_SPIKE_SPIKE2DATA, use
%   [freq] = spiketriggeredspectrum(cfg, data)
% When spikes from SPIKE struct need to be directly combined with continuous representation, use as
%   [freq] = ft_spike_triggeredspectrum(cfg,data,spike)
% Important: in that case data.time and spike.trialtime should be matched!
%
%   Configurations:
%   cfg.timwin       = [begin end], time around each spike (default = [-0.1 0.1])
%   cfg.foilim       = [begin end], frequency band of interest (default = [0 150])
%   cfg.taper        = 'hanning' or many others, see WINDOW (default = 'hanning'). 
%   cfg.tapsmofrq    = number, the amount of spectral smoothing through
%                      multi-tapering. Note that 4 Hz smoothing means
%                      plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
%   cfg.taperopt     = additional inputs going to WINDOW func as cell array, e.g.
%                      cfg.taperopt = {5, -35} with cfg.taper = 'taylorwin'.
%                      Note: multitapering rotates phases (no problem for consistency)
%   cfg.spikechannel = string, name of single spike channel to trigger on. See FT_CHANNELSELECTION
%   cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                      see FT_CHANNELSELECTION for details
%   cfg.feedback     = 'no', 'text', 'textbar', 'gui' (default = 'no')
%   cfg.algorithm    = 'vectorized' (default) or 'loop'. If the window around each spike
%                      is short, a moderate a speed up of 5-10 times is typical using vectorized 
%                      code. However, for long windows (e.g. 0.5-1 second), vectorized
%                      algorithm is typically 50-100 times faster because it calls FFT only 
%                      once per trial (instead of nSpikes times).
%                      In case an lfp segment contains NaNs, we use specest_nanfft.
%   cfg.borderspikes = 'yes' (default) or 'no'. If 'yes', we process the spikes falling at the
%                      border using the same taper for the LFP that falls within the trial.
%
% If the triggered spike leads a spike in another channel, then the angle
% of the Fourier spectrum of that other channel will be negative. NOTE that
% this should be checked for consistency.
% % $Id$

% Copyright (C) 2008-2011, Robert Oostenveld, Martin Vinck

% set the defaults
if ~isfield(cfg, 'timwin'),         cfg.timwin         = [-0.1 0.1];   end
if ~isfield(cfg, 'foilim'),         cfg.foilim         = [0 150];      end
if ~isfield(cfg, 'taper'),          cfg.taper          = 'hanning';    end
if ~isfield(cfg, 'channel'),        cfg.channel        = 'all';        end
if ~isfield(cfg, 'spikechannel'),   cfg.spikechannel   = [];           end
if ~isfield(cfg, 'feedback'),       cfg.feedback       = 'no';         end
if ~isfield(cfg, 'algorithm'),      cfg.algorithm      = 'loop';       end
if ~isfield(cfg, 'taperopt'),       cfg.taperopt       = [];           end
if ~isfield(cfg,'borderspikes'),    cfg.borderspikes   = 'yes';        end
  
if strcmp(cfg.taper, 'sine')
  error('sorry, sine taper is not yet implemented');
end

% determine which algorithm to use
if strcmp(cfg.algorithm,'loop')
  doLoop = 1;
elseif strcmp(cfg.algorithm,'vectorized')
  doLoop = 0;
else
  error('unrecognized option for cfg.algorithm')
end

doBorderSpikes = strcmp(cfg.borderspikes,'yes');

% check whether a third input is specified
hasSpike = nargin==3;

% do the selection of the spikechannels
if hasSpike,
  spikelabel = spike.label;
else
  spikelabel = data.label;
end
  
% determine the channels to be averaged
cfg.channel = ft_channelselection(cfg.channel, data.label);
chansel     = match_str(data.label, cfg.channel); % selected channels
nchansel    = length(cfg.channel);                % number of channels

% get the spikechannels
cfg.spikechannel = ft_channelselection(cfg.spikechannel, spikelabel);
spikesel         = match_str(spikelabel, cfg.spikechannel);
nspikesel        = length(cfg.spikechannel); % number of spike channels
if nspikesel~=1, error('MATLAB:ft_spike_triggeredspectrum:cfg:spiketriggeredspectrum:notOneSelected',...
                    'Onle one spike channel can be selected at a time'); 
end

nTrials    = length(data.trial); % number of trials

% calculate number of samples left and right from spike and construct the taper
begpad = round(cfg.timwin(1)*data.fsample); % number of samples before spike
endpad = round(cfg.timwin(2)*data.fsample); % number of samples after spike
numsmp = endpad - begpad + 1;               % number of samples of segment
if ~strcmp(cfg.taper,'dpss')
  if iscell(cfg.taperopt)
    taper  = window(cfg.taper, numsmp, cfg.taperopt{:});
  else
    taper  = window(cfg.taper, numsmp);
  end
  taper  = taper./norm(taper);
else
  % not implemented yet: keep tapers, or selecting only a subset of them.
  taper  = dpss(numsmp, cfg.tapsmofrq);
  taper  = taper(:,1:end-1);            % we get 2*NW-1 tapers
  taper  = sum(taper,2)./size(taper,2); % using the linearity of multitapering
end

% frequency selection
freqaxis = linspace(0, data.fsample, numsmp);
fbeg     = nearest(freqaxis, cfg.foilim(1));
fend     = nearest(freqaxis, cfg.foilim(2));
freqs    = freqaxis(fbeg:fend)'; % column vector that we use for rephasing later
nFreqs   = fend - fbeg + 1;
cfg.foilim(1) = freqaxis(fbeg); % update the configuration to account for rounding off differences
cfg.foilim(2) = freqaxis(fend);

% make a representation of the spike, this is used for the phase rotation
spk = zeros(numsmp,1);
spk(1-begpad) = 1;
spike_fft = fft(spk);
spike_fft = spike_fft(fbeg:fend);
spike_fft = spike_fft./abs(spike_fft);
rephase   = conj(spike_fft);

% preallocate the outputs
spectrum    = cell(1,nTrials);               
spiketime   = cell(1,nTrials);
spiketrial  = cell(1,nTrials);

% compute the spectra
for iTrial = 1:nTrials

    if hasSpike 
        % get all the samples at once without using loops
        hasTrial = spike.trial{spikesel} == iTrial;
        ts       = spike.time{spikesel}(hasTrial);
        spikesmp = nearest_nd(data.time{iTrial},ts); % this gives unique samples back     

        % calculate the phase shift based on difference between 'exact' and sample time
        ts       = ts(:);
        smptime  = data.time{iTrial}(spikesmp);
        shift    = ts - smptime(:);

    elseif ~hasSpike
        spikesmp = find(data.trial{iTrial}(spikesel,:)); % simply find the LFP segment belonging to the time-stamp here.
        spikecnt = data.trial{iTrial}(spikesel,spikesmp);

        % instead of doing the bookkeeping of double spikes below, replicate the double spikes by looking at spikecnt
        sel = find(spikecnt>1);
        tmp = zeros(1,sum(spikecnt(sel)));
        n   = 1;
        for j=1:length(sel)
            for k=1:spikecnt(sel(j))
                tmp(n) = spikesmp(sel(j));
                n = n + 1;
            end
        end

        % remove the double spikes
        spikesmp(sel) = [];                     % remove the double spikes
        spikesmp = sort([spikesmp tmp]);        % add the double spikes as replicated single spikes
    end

    % determine the beginning and end of the segments
    begsmp = spikesmp + begpad;
    endsmp = spikesmp + endpad;

    % determine which spikes do not have their windows at the edges of the data
    nSpikes    = length(begsmp);
    vldSpikes    = begsmp>=1&endsmp<=size(data.trial{iTrial},2); % determine the valid spikes
    earlySpikes  = begsmp<1 & spikesmp>=1;
    lateSpikes   = spikesmp<=size(data.trial{iTrial},2) & endsmp>size(data.trial{iTrial},2);
    if doBorderSpikes
      vldSpikes  = find(vldSpikes(:) | earlySpikes(:) | lateSpikes(:));
      begsmp(earlySpikes) = 1;
      endsmp(earlySpikes) = numsmp;
      endsmp(lateSpikes)  = size(data.trial{iTrial},2);
      begsmp(lateSpikes)  = size(data.trial{iTrial},2)-numsmp+1;
    else
      vldSpikes  = find(vldSpikes);
    end
      
    nVld       = length(vldSpikes);

    % we want to process the incomplete spikes, also with hann window
    
    % determine the new begsmp and endsmp for the spikes that fall around the border, we should not
    % discard these, esp. not for low frequencies.
    fprintf('processing trial %d of %d (%d spikes)\n', iTrial, nTrials, sum(nSpikes));

    % store the spiketimes and spiketrials, constructing a type of SPIKE format
    if ~hasSpike % otherwise we can immediate copy to FREQ.origtime and FREQ.origtrial
        spiketime{iTrial}  = data.time{iTrial}(spikesmp);
        spiketrial{iTrial} = iTrial*ones(1,nSpikes);
    end
    
    % preallocate the spectrum that we create, also the spikes that will have NaNs
    spectrum{iTrial} = NaN(nSpikes, nchansel, nFreqs); 
    rephaseMat = rephase(:,ones(1,nVld),ones(1,nchansel));
    if doBorderSpikes && sum(earlySpikes)>0
      for k = find(earlySpikes)        
        spk = zeros(numsmp,1);
        spk(spikesmp(k)) = 1;
        spike_fft = fft(spk);
        spike_fft = spike_fft(fbeg:fend);
        spike_fft = spike_fft./abs(spike_fft);
        reph   = conj(spike_fft);
        rephaseMat(:,k,:) = repmat(reph,[1 1 nchansel]);
      end
    end
    if doBorderSpikes && sum(lateSpikes)>0
      n  = size(data.trial{iTrial},2);
      bg = n-numsmp; 
      for k = find(lateSpikes)        
        spk = zeros(numsmp,1);
        spk(spikesmp(k)-bg) = 1;
        spike_fft = fft(spk);
        spike_fft = spike_fft(fbeg:fend);
        spike_fft = spike_fft./abs(spike_fft);
        reph   = conj(spike_fft);
        rephaseMat(:,k,:) = repmat(reph,[1 1 nchansel]);
      end
    end
    
    
    if nVld<1,continue,end % continue to the next trial if this one does not contain valid spikes

    if ~doLoop
      % get the lfp, and reshape it such that we can index it at once.
      lfp     = data.trial{iTrial}(chansel,:);
      lfp     = reshape(lfp', [size(lfp,2) 1 nchansel]); % reshape to by time-by-1-by-nchansel

      % form the indices to index the LFP at once
      segment = linspacevec(begsmp(vldSpikes),endsmp(vldSpikes),numsmp)';
      
      % cut out all windows at once to form a large 2-D matrix
      segment = lfp(segment,:,:);
      segment = reshape(segment,[numsmp nVld nchansel]);
      segmentMean = repmat(nanmean(segment,1),[numsmp 1 1]); % nChan x Numsmp
      segment     = segment - segmentMean; % LFP has average of zero now (no DC)
      segmentMean = []; % remove to save memory

      lfp     = [];

      hasNan = find(any(isnan(segment))); % this will be indx for nSpikes-by-nChans
      [row,col] = ind2sub([nVld nchansel], hasNan);

      % replace the parts in segment where there were nans in the data
      %my_fft = @specest_nanfft;
      my_fft = @specest_nanfft_;
      nNans  = length(row);
      fftNan = zeros(numsmp,nNans);
      time  = randn(size(segment)); % this is actually not used (MV: why is it there then?)
      for iNan = 1:nNans
        fftNan(:,iNan) = my_fft(segment(:,row(iNan),col(iNan)).*taper, time);
      end        

      % fourier transform the data matrix 
      segment = fft(segment .* taper(:,ones(1,nVld), ones(1,nchansel)));

      % replace the part of the segment with NaNs now
      for iNan = 1:nNans
        segment(:,row(iNan),col(iNan)) = fftNan(:,iNan);
      end

      % select the desired output frequencies and normalize
      segment = segment(fbeg:fend,:,:) ./ sqrt(numsmp/2);

      % rotate the phase at each frequency to correct for segment t=0 not being at the first sample
      if hasSpike
          % phase rotation according to difference data timeaxis and raw timestamps
          shift = shift(vldSpikes)';  
          shiftrephase = shift(ones(nFreqs,1),:,ones(1,nchansel)).*freqs(:,ones(1,nVld),ones(1,nchansel));
          shiftrephase = exp(i*2*pi*shiftrephase);

          % rephase the fourier transform
          segment = segment .* (rephase(:,ones(1,nVld),ones(1,nchansel)) .* shiftrephase); 
      else
          segment = segment .* (rephase(:,ones(1,nVld),ones(1,nchansel)));
      end

      % collect the results
      spectrum{iTrial}(vldSpikes,1:nchansel,1:nFreqs) = permute(segment,[2 3 1]);
    else        
      my_fft = @specest_nanfft;
      ft_progress('init', cfg.feedback, 'spectrally decomposing data around spikes');
      reph = repmat(rephase,1,nchansel).';
      for iSpike = vldSpikes
          ft_progress(iTrial/nTrials, 'spectrally decomposing data around spike %d of %d\n', iSpike, length(spikesmp));
          begsmp = spikesmp(iSpike) + begpad;
          endsmp = spikesmp(iSpike) + endpad;
          segment = data.trial{iTrial}(chansel,begsmp:endsmp);

          % substract the DC component from every segment, to avoid any leakage of the taper
          segmentMean = repmat(nanmean(segment,2),1,numsmp); % nChan x Numsmp
          segment     = segment - segmentMean; % LFP has average of zero now (no DC)

          % taper the data segment around the spike and compute the fft
          time  = randn(size(segment)); % this is actually not used (MV: why is it there then?)                                
          segment_fft = my_fft(segment * sparse(diag(taper)),time);

          % select the desired output frequencies and normalize
          segment_fft = segment_fft(:,fbeg:fend) ./ sqrt(numsmp/2);

          % rotate the phase at each frequency to correct for segment t=0 not being at the first sample
          if hasSpike
              % phase rotation according to difference data timeaxis and raw timestamps
              spikerephase = exp(2*pi*1i*shift(iSpike)*freqs); % use fft shift property
              segment_fft = segment_fft.* reph .* spikerephase(:,ones(1,nchansel)).';
          else
              segment_fft = segment_fft.* reph;
          end

          % store the result for this spike in this trial
          spectrum{iTrial}(iSpike,:,:) = segment_fft;

      end % for each spike in this trial
      ft_progress('close');
    end
end % for each trial

% collect the results
freq.label          = data.label(chansel);
freq.freq           = freqaxis(fbeg:fend);
freq.dimord         = 'rpt_chan_freq';
freq.fourierspctrm  = cat(1, spectrum{:});
if ~hasSpike
  freq.time     = cat(2,spiketime{:})';  
  freq.trial    = cat(2,spiketrial{:})'; 
  for iTrial = 1:nTrials
    freq.trialtime(iTrial,:) = [min(data.time{iTrial}) max(data.time{iTrial})];
  end
else
  freq.time       = spike.time{spikesel}(:);  
  freq.trial      = spike.trial{spikesel}(:); 
  freq.trialtime  = spike.trialtime;
end  

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, ind] = dbstack;
  cfg.version.name = st(ind);
end
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
freq.cfg = cfg;

function [indx] = nearest_nd(x,y)

% NEAREST return the index of an n-d matrix to an n-d matrix.
% 
% [indx] = nearest_nd(x, y)
%
% Inputs: 
%   X can be a n-d matrix of any size (scalar, vector, n-d matrix).
%   Y can be an n-d matrix of any size (scalar, vector, n-d matrix).
%   
% If Y is larger than any X, we return the last index that the maximum value of X occurred.
% Otherwise, we return the first occurence of the nearest X.
%
% If Y contains NaNs, we return a NaN for every NaN in Y.
%
% Outputs:
%   INDX is a vector of size Y and contains the indices of the values in X that are
%   closest to the respective value of Y. INDX is a linear index, such that x(INDX) gives
%   the nearest values of X to Y. To convert INDX to subscripts, see IND2SUB.
%
% Example:
% y = 100*rand(5,10);
% x = 1:100;
% indx = nearest_nd(x,y);
% disp(max(max(y-x(indx))))
%
% x = reshape([1:100 100*ones(1,20)],[20 6]);
% y = rand(5,10)*max(x(:));
% y(1) = 100.2;
% indx = nearest_nd(x,y);
% disp(max(max(y-x(indx))))
% Note that indx(1) returns the last occurence and indx is a linear index to 2-d X.
% 
% Demonstrate the use with a 3-D x and y array.
% x = reshape([1:100],[10 5 2]);
% y = rand(5,10,2)*max(x(:));
% indx = nearest_nd(x,y);
% d = x(indx)-y; disp(max(d(:)));
% y(1) = NaN;
% indx = nearest_nd(x,y);
% disp(indx(1))
%

% Copyright, Martin Vinck, 2009.

% store the sizes of x and y, this is used to reshape INDX later on
szY = size(y);

% vectorize both x and y
x = x(:);
y = y(:);

% from now on we can treat X and Y as vectors
nY = length(y);
nX = length(x);
hasNan = isnan(y); % indices with nans in y, indx(hasNaN) will be set to NaN later.

if nX==1, 
  indx = ones(1,nY);  % only one x value, so nearest is always only element
else
  if nY==1 % in this case only one y value, so use the old NEAREST code
    if y>max(x)
      % return the last occurence of the nearest number
      [dum, indx] = max(flipud(x));
      indx = nX + 1 - indx;
    else
      % return the first occurence of the nearest number
      [mindist, indx] = min(abs(x(:) - y));
    end
  else
    if any(y>max(x))
      % for these return the last occurence of every number as in NEAREST
      indx = zeros(1,nY);
      i = y>max(x);
      [dum,indx] = max(flipud(x));
      indx(i)       = nX + 1 - indx;    
      % for the rest return the first occurence of every number
      x = x(:); 
      y = y(~i)';
      xRep = x(:,ones(1,length(y)));
      yRep = y(ones(nX,1),:);
      [mindist,indx(~i)] = min(abs(xRep-yRep));
    else
      x = x(:); 
      y = y';
      xRep = x(:,ones(1,nY));
      yRep = y(ones(nX,1),:);
      [mindist,indx] = min(abs(xRep-yRep));
    end        
  end
end
% return a NaN in INDX for a NaN in Y
indx(hasNan) = NaN;

% reshape the indx back to the y format
if (sum(szY>1)>1 || length(szY)>2) % in this case we are dealing with a matrix
  indx = reshape(indx,[szY]);
end

function y = linspacevec(d1, d2, n)

%LINSPACEVEC Linearly spaced matrix with different stop and end values.
% Example:
% d1 = 1:10;
% d2   = 11:2:30;
% n = 10;
% indx = linspacevec(d1,d2,n)

%   Copyright 2008 Vinck

if nargin == 2
    n = 100;
end

d1 = d1(:);
d2 = d2(:);

d = d2-d1;
if numel(d1)==1
    n = double(n);
    y = [d1+(0:n-2)*d/(floor(n)-1) d2];
else    
    D = d(:,ones(n-1,1));
    v = (0:n-2)';
    N = v(:,ones(length(d1),1))';
    D1 = d1(:,ones(n-1,1));
    y = [D1+N.*D/(n-1) d2];
end

  