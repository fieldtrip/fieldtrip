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
%   cfg.taper        = 'hanning' or many others, see WINDOW (default = 'hanning'). Multi-tapering is
%                       not implemented since this will rotate the phases. 
%   cfg.spikechannel = string, name of single spike channel to trigger on. See FT_CHANNELSELECTION
%   cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                      see FT_CHANNELSELECTION for details
%   cfg.feedback     = 'no', 'text', 'textbar', 'gui' (default = 'no')
%
% If the triggered spike leads a spike in another channel, then the angle
% of the Fourier spectrum of that other channel will be negative. NOTE that
% this should be checked for consistency.
%

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: spiketriggeredspectrum.m,v $
% Revision Martin Vinck Oct 2010:
% Input can now also be spike, because this will be more accurate for high frequencies.
% Changed output to be compatible with SPIKE type representation, adding field trialtime
% Removed checking of spike channel inside function; this should be solved with cfg at least,
% because sometimes we do not use spikes but events for example. Also depends a lot on recording
% system
%
% Revision 1.8:
% - Allowing for any window (hanning was entered always).
% - Keeping the NaNs, because this makes it much easier to link STS structures
%   with different frequency content.
% - New spike channel autodetection function, and only checking for selected spike channel.
% - Added an optional vectorized algorithm, for longer windows around the spikes this gives
%   a great speed-up (up to factor 100) without loss of numerical precision.
% - Allowed for SPIKE input, and a rephasing of the spike phase by using the exact spike time
%   (which is slightly more accurate than using the sample time, can be relevant for high freqs).
%
% Revision 1.7  2008/10/01 08:23:05  roboos
% use specest_nanfft as default
%
% Revision 1.6  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.5  2008/07/15 18:17:57  roboos
% fixed documentation
%
% Revision 1.4  2008/04/09 14:39:22  roboos
% remove trials without data (too close to trial border) from the output
%
% Revision 1.3  2008/04/09 14:23:15  roboos
% replicate double spikes, for keeptrials
% give error if spike count >5
%
% Revision 1.2  2008/03/18 21:58:18  roboos
% added normalization of the spectrum by the number of samples, and rotation of the phase to correct for t=0 being in the middle of the segment instead of at the begin
%
% Revision 1.1  2008/03/17 16:42:35  roboos
% initial implementation, phase rotation still needs to be built in
%
% add cell evaluation for speed up

% set the defaults
if ~isfield(cfg, 'timwin'),         cfg.timwin         = [-0.1 0.1];   end
if ~isfield(cfg, 'foilim'),         cfg.foilim         = [0 150];      end
if ~isfield(cfg, 'taper'),          cfg.taper          = 'hanning';    end
if ~isfield(cfg, 'channel'),        cfg.channel        = 'all';        end
if ~isfield(cfg, 'spikechannel'),   cfg.spikechannel   = [];           end
if ~isfield(cfg, 'feedback'),       cfg.feedback       = 'no';         end

if strcmp(cfg.taper, 'dpss') || strcmp(cfg.taper, 'sine') 
    error('sorry, multitapering is not yet implemented');
end

% check whether a third input is specified
hasSpike = nargin==3;

% do the selection of the spikechannels
if hasSpike,
  label = spike.label;
else
  label = data.label;
end
  
% determine the channels to be averaged
cfg.channel = ft_channelselection(cfg.channel, label);
chansel     = match_str(label, cfg.channel); % selected channels
nchansel    = length(cfg.channel);                % number of channels

% get the spikechannels
cfg.spikechannel = ft_channelselection(cfg.spikechannel, label);
spikesel         = match_str(label, cfg.spikechannel);
nspikesel        = length(cfg.spikechannel); % number of spike channels
if nspikesel~=1, error('MATLAB:ft_spike_triggeredspectrum:cfg:spiketriggeredspectrum:notOneSelected',...
                    'Onle one spike channel can be selected at a time'); 
end

nTrials    = length(data.trial); % number of trials

% calculate number of samples left and right from spike and construct the taper
begpad = round(cfg.timwin(1)*data.fsample); % number of samples before spike
endpad = round(cfg.timwin(2)*data.fsample); % number of samples after spike
numsmp = endpad - begpad + 1;               % number of samples of segment
taper  = window(cfg.taper, numsmp);         % create a hanning taper
taper  = taper./norm(taper);                % normalize the taper

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
rephase   = rephase(:,ones(1,nchansel)).'; % bring it in the shape we will need later on

% preallocate the outputs
spectrum    = cell(1,nTrials);               
spiketime   = cell(1,nTrials);
spiketrial  = cell(1,nTrials);

% compute the spectra
for iTrial = 1:nTrials

    if hasSpike 
        % get all the samples at once without using loops
        hasTrial = spike.trial{spikesel} == iTrial;
        ts       = spike.timestamp{spikesel}(hasTrial);
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
    vldSpikes  = find(begsmp>=1&endsmp<=size(data.trial{iTrial},2)); % determine the valid spikes
    nVld       = length(vldSpikes);

    fprintf('processing trial %d of %d (%d spikes)\n', iTrial, nTrials, sum(nSpikes));

    % store the spiketimes and spiketrials, constructing a type of SPIKE format
    if ~hasSpike % otherwise we can immediate copy to FREQ.origtime and FREQ.origtrial
        spiketime{iTrial}  = data.time{iTrial}(spikesmp);
        spiketrial{iTrial} = iTrial*ones(1,nSpikes);
    end
    
    % preallocate the spectrum that we create, also the spikes that will have NaNs
    spectrum{iTrial} = NaN(nSpikes, nchansel, nFreqs); 

    if nVld<1,continue,end % continue to the next trial if this one does not contain valid spikes

    my_fft = @specest_nanfft;
    progress('init', cfg.feedback, 'spectrally decomposing data around spikes');
    for iSpike = vldSpikes
        progress(iTrial/nTrials, 'spectrally decomposing data around spike %d of %d\n', iSpike, length(spikesmp));
        begsmp = spikesmp(iSpike) + begpad;
        endsmp = spikesmp(iSpike) + endpad;
        segment = data.trial{iTrial}(chansel,begsmp:endsmp);

        % taper the data segment around the spike and compute the fft
        segment_fft = my_fft(segment * sparse(diag(taper)));

        % select the desired output frequencies and normalize
        segment_fft = segment_fft(:,fbeg:fend) ./ sqrt(numsmp/2);

        % rotate the phase at each frequency to correct for segment t=0 not being at the first sample
        if hasSpike
            % phase rotation according to difference data timeaxis and raw timestamps
            spikerephase = exp(2*pi*1i*shift(iSpike)*freqs); % use fft shift property
            segment_fft = segment_fft.* rephase .* spikerephase(:,ones(1,nchansel)).';
        else
            segment_fft = segment_fft.* rephase;
        end

        % store the result for this spike in this trial
        spectrum{iTrial}(iSpike,:,:) = segment_fft;

    end % for each spike in this trial
    progress('close');
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
    freq.trialtime(iTrial,:) = [min(data.trial{iTrial}) max(data.trial{iTrial})];
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

  