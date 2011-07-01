function [spectrum, freqoi, timeoi] = ft_specest_convol(dat, time, varargin)

% SPECEST_CONVOL performs wavelet convolution in the time domain by
% convolution with Morlet's wavelets.
%
% Use as
%   [spectrum,freqoi,timeoi] = specest_convol(dat,time,...)
% where
%   dat      = matrix of chan*sample
%   time     = vector, containing time in seconds for each sample
%   spectrum = matrix of chan*freqoi*timeoi of fourier coefficients
%   freqoi   = vector of frequencies in spectrum
%   timeoi   = vector of timebins in spectrum
%
% Optional arguments should be specified in key-value pairs and can include:
%   timeoi        = vector, containing time points of interest (in seconds, analysis window will be centered around these time points)
%   freqoi        = vector, containing frequencies (in Hz)
%   waveletwidth  = number, 'width' of wavelets expressed in cycles (default = 7)
%
%  OPTION DOWNSAMPLE: this looks like something that would fit in better in the (freqanalysis) wrapper I think
%  OPTION LATENCY, WHAT TO DO WITH THIS?
%  HOW TO MAKE CONSISTENT freqoi'S OVER TRIALS (e.g. multiple runs of this function)? ZERO-PADDING NOT AN OPTION... YET WE SHOULD STILL ONLY HAVE freqoi'S THAT ACTUALLY MATCH THE FREQUENCIES IN OUTPUT
%
% See also SPECEST_MTMFFT, SPECEST_MTMCONVOL, SPECEST_HILBERT, SPECEST_NANFFT, SPECEST_MVAR, SPECEST_WAVELET

% Copyright (C) 2010, Donders Institute for Brain, Cognition and Behaviour
%
% $Log$

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'waveletwidth','pad','timeoi','freqoi','polyremoval'});
timeoi        = keyval('timeoi',        varargin); if isempty(timeoi),       timeoi       = 'all';    end
freqoi        = keyval('freqoi',        varargin); if isempty(freqoi),       freqoi       = 'all';    end
waveletwidth  = keyval('waveletwidth',  varargin); if isempty(waveletwidth), waveletwidth = 7;        end
polyorder     = keyval('polyremoval', varargin); if isempty(polyorder), polyorder = 1; end

% Set n's
[nchan,ndatsample] = size(dat);

% Remove polynomial fit from the data -> default is demeaning
if polyorder >= 0
  dat = ft_preproc_polyremoval(dat, polyorder, 1, ndatsample);
end

% Determine fsample and set total time-length of data
fsample = 1/(time(2)-time(1));
dattime = ndatsample / fsample; % total time in seconds of input data
endnsample = ndatsample;  % for consistency with mtmconvol and mtmfft
endtime    = dattime;     % for consistency with mtmconvol and mtmfft

% Set freqboi and freqoi
if isnumeric(freqoi) % if input is a vector
  freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1;
  freqboi   = unique(freqboi);
  freqoi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output
elseif strcmp(freqoi,'all')
  freqboilim = round([0 fsample/2] ./ (fsample ./ endnsample)) + 1;
  freqboi    = freqboilim(1):1:freqboilim(2);
  freqoi     = (freqboi-1) ./ endtime;
end
nfreqboi   = length(freqboi);
nfreqoi = length(freqoi);

% Set timeboi and timeoi
offset = round(time(1)*fsample);
if isnumeric(timeoi) % if input is a vector
  timeboi  = round(timeoi .* fsample - offset) + 1;
  ntimeboi = length(timeboi);
  timeoi   = round(timeoi .* fsample) ./ fsample;
elseif strcmp(timeoi,'all') % if input was 'all'
  timeboi  = 1:length(time);
  ntimeboi = length(timeboi);
  timeoi   = time;
end

% compute wavelet family
wavfam = waveletfam(freqoi,fsample,waveletwidth);

% compute spectrum by convolving the wavelets with the data
spectrum = zeros(nchan, nfreqoi, nsample);
for ifreqoi = 1:nfreqoi
  wavelet = wavfam{ifreqoi};
  for ichan = 1:nchan
    spectrum(ichan,ifreqoi,:) = conv(dat(ichan,:),  wavelet, 'same');
  end
  % pad the edges with nans to indicate that the wavelet was not fully immersed in the data THERE ARE NO NANS ADDED IN FREQANALYSIS_TFR?
  nanpad = ceil(length(wavfam{ifreqoi})/2);
  % the padding should not be longer than the actual data
  nanpad = min(nanpad, nsample);
  begnanpad = 1:nanpad;
  endnanpad = (nsample-nanpad+1):nsample;
  spectrum(:,ifreqoi,begnanpad) = nan;
  spectrum(:,ifreqoi,endnanpad) = nan;
end

% select the samples for the output
spectrum = spectrum(:,:,timeboi);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for waveletanalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = waveletfam(foi,fsample,waveletwidth)
dt = 1/fsample;
for k=1:length(foi)
  sf  = foi(k)/waveletwidth;
  st  = 1/(2*pi*sf);
  toi = -3.5*st:dt:3.5*st;
  A   = 1/sqrt(st*sqrt(pi));
  M{k}= A*exp(-toi.^2/(2*st^2)).*exp(i*2*pi*foi(k).*toi);
end



% state         = keyval('state',         varargin);
% waveletwidth  = keyval('waveletwidth',  varargin); if isempty(waveletwidth), waveletwidth = 7; end
% downsample    = keyval('downsample',    varargin);
% time          = keyval('time',          varargin);
% toi           = keyval('toi',           varargin);
%
% % FIXME add width
%
% if ~isempty(downsample) && ~isempty(toi)
%   error('the downsample and toi options are mutually exclusive');
% end
%
% nchans   = size(dat,1);
% nsamples = size(dat,2);
%
% if isempty(time)
%   if ~isempty(toi)
%     error('specification of toi without time is not permitted');
%   end
%   time = (0:(nsamples-1))/fsample;
% end
%
% if ~isempty(state) && isequal(state.fsample, fsample)  && isequal(state.foi, foi) && isequal(state.waveletwidth, waveletwidth)
%   % reuse the wavelet family from the previous state
%   disp('using previous state');
%   M = state.M;
% else
%   if ~isempty(state)
%     state = [];
%     warning('recomputing the state');
%   end
%   % recompute the wavelet family
%   M = waveletfam(foi,fsample,waveletwidth);
% end
%
% % compute freq by convolving the wavelets with the data
% nfreq = length(foi);
% spectrum = zeros(nchans, nsamples, nfreq);
%
% for f=1:nfreq
%   wavelet = M{f};
%   for c=1:nchans
%     spectrum(c,:,f) = conv(dat(c,:),  wavelet, 'same');
%   end
%   % pad the edges with nans to indicate that the wavelet was not fully immersed in the data
%   pad = ceil(length(M{f})/2);
%   % the padding should not be longer than the actual data
%   pad = min(pad, nsamples);
%   begpad = 1:pad;
%   endpad = (nsamples-pad+1):nsamples;
%   spectrum(:,begpad,f) = nan;
%   spectrum(:,endpad,f) = nan;
% end
%
% if ~isempty(downsample)
%   % this would be done in case of downsample
%   tbin = 1:downsample:nsamples;
% elseif ~isempty(toi)
%   % this would be done in case of toi/time specification
%   tbin = zeros(size(toi));
%   for i=1:length(tbin)
%     tbin = nearest(time, toi(i));
%   end
% else
%   tbin = 1:nsamples;
% end
%
% % select the samples for the output
% spectrum = spectrum(:,tbin,:);
%
% % remember the wavelets so that they can be reused in a subsequent call
% state.fsample = fsample;
% state.foi     = foi;
% state.waveletwidth = waveletwidth;
% state.M       = M;
%
% % %convolves data in chan by time matrix (from one trial) with 1 wavelet for
% % %each specified frequency of interest
% % spectrum = zeros(size(data,1),length(foi), size(data,2));
% % for k=1:size(data,1) %nchans
% %   for j=1:length(foi)
% %     cTmp = conv(data(k,:),M{j});
% %     cTmp = 2*(abs(cTmp).^2)/fsample;
% %     cTmp = cTmp(ceil(length(M{j})/2):length(cTmp)-floor(length(M{j})/2));
% %     spectrum(k,j,:) = cTmp;
% %   end
% %
% % end
%
% %output should be chan by freq by time
