function [spectrum, foi, toi] = specest_tfr(dat, time, varargin)

% SPECEST_TFR performs wavelet convolution in the time domain by convolution with Morlet's wavelets.
%
%
% Use as
%   [spectrum,foi,toi] = specest_tfr(dat,time,...)  
%
%   dat      = matrix of chan*sample
%   time     = vector, containing time in seconds for each sample
%   spectrum = matrix of chan*foi*toi of fourier coefficients
%   foi      = vector of frequencies in spectrum
%   toi      = vector of timebins in spectrum
%
%
%
%
% Optional arguments should be specified in key-value pairs and can include:
%   toi           = vector, containing time points of interest (in seconds, analysis window will be centered around these time points)
%   foi           = vector, containing frequencies (in Hz)
%   waveletwidth  = number, 'width' of wavelets expressed in cycles (default = 7)   
%
%
%
%  OPTION DOWNSAMPLE: this looks like something that would fit in better in the wrapper I think
%
%  OPTION LATENCY, WHAT TO DO WITH THIS? 
%
%  HOW TO MAKE CONSISTENT FOI'S OVER TRIALS? ZERO-PADDING NOT AN OPTION... YET WE SHOULD STILL ONLY HAVE FOI'S THAT ACTUALLY MATCH THE FREQUENCIES IN OUTPUT
%

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'waveletwidth','pad','toi','foi'});
toi           = keyval('toi',           varargin); if isempty(toi),          toi          = 'max';    end
foi           = keyval('foi',           varargin); if isempty(foi),          foi          = 'max';    end
waveletwidth  = keyval('waveletwidth',  varargin); if isempty(waveletwidth), waveletwidth = 7;        end


% Set n's
[nchan,nsample] = size(dat);


% Determine fsample
fsample = nsample / (time(end) - time(1));


% Set fboi and foi and default tapsmofrq
pad = (time(end)-time(1)); % PADDING DOES NOT EXIST IN TFR, PERHAPS MAKE A NEW VAR-NAME TO BE USED IN ALL SPECEST, LIKE DATLENGHT OR SOMETHING, FOR CONSISTENCY
if isnumeric(foi) % if input is a vector
  fboi   = round(foi ./ (fsample ./ (pad * fsample))) + 1;
  nfboi  = length(fboi);
  foi    = (fboi-1) ./ pad; % boi - 1 because 0 Hz is included in fourier output..... is this going correctly?
elseif strcmp(foi,'max') % if input was 'max'
  fboilim = round([0 fsample/2] ./ (fsample ./ (pad * fsample))) + 1;
  fboi    = fboilim(1):1:fboilim(2);
  nfboi   = length(fboi);
  foi     = (fboi-1) ./ pad;
end
% check for foi = 0 and remove it, there is no wavelet for foi = 0
if any(foi==0)
  foi(foi==0) = [];
end
nfoi = length(foi);


% Set tboi and toi
if isnumeric(toi) % if input is a vector
  tboi  = round(toi .* fsample) + 1;
  ntboi = length(tboi);
  toi   = round(toi .* fsample) ./ fsample;
elseif strcmp(toi,'max') % if input was 'max'
  tboi  = 1:length(time);
  ntboi = length(tboi);
  toi   = time;
end



% compute wavelet family
wavfam = waveletfam(foi,fsample,waveletwidth);


% compute spectrum by convolving the wavelets with the data
spectrum = zeros(nchan, nfoi, nsample);
for ifoi = 1:nfoi
  wavelet = wavfam{ifoi};
  for ichan = 1:nchan
    spectrum(ichan,ifoi,:) = conv(dat(ichan,:),  wavelet, 'same');
  end
  % pad the edges with nans to indicate that the wavelet was not fully immersed in the data THERE ARE NO NANS ADDED IN FREQANALYSIS_TFR?
  nanpad = ceil(length(wavfam{ifoi})/2);
  % the padding should not be longer than the actual data
  nanpad = min(nanpad, nsample);
  begnanpad = 1:nanpad;
  endnanpad = (nsample-nanpad+1):nsample;
  spectrum(:,ifoi,begnanpad) = nan;
  spectrum(:,ifoi,endnanpad) = nan;
end

% select the samples for the output
spectrum = spectrum(:,:,tboi);





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













