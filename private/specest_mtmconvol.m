function [spectrum,ntaper,freqoi,timeoi] = specest_mtmconvol(dat, time, varargin)

% SPECEST_MTMCONVOL performs wavelet convolution in the time domain by multiplication in the frequency domain
%
%
% Use as
%   [spectrum,ntaper,freqoi,timeoi] = specest_mtmconvol(dat,time,...) 
%
%   dat      = matrix of chan*sample
%   time     = vector, containing time in seconds for each sample
%   spectrum = matrix of taper*chan*freqoi*timeoi of fourier coefficients
%   ntaper   = vector containing number of tapers per element of freqoi
%   freqoi   = vector of frequencies in spectrum
%   timeoi   = vector of timebins in spectrum
%
%
%
%
% Optional arguments should be specified in key-value pairs and can include:
%   taper     = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%   pad       = number, indicating time-length of data to be padded out to in seconds         
%   timeoi    = vector, containing time points of interest (in seconds)
%   timwin    = vector, containing length of time windows (in seconds)
%   freqoi    = vector, containing frequencies (in Hz)
%   tapsmofrq = vector, the amount of spectral smoothing through multi-tapering. Note: 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box
%
%
%
%
%
% FFT SPEED NOT YET OPTIMIZED (e.g. matlab version, transpose or not)
% IF FREQOI CONTAINS 0 (we should either remove it or allow the creation of a wavelet which is a straigth line, it's removed now)
% SHOULD FREQOI = 'ALL' BE REMOVED OR NOT?
%
% See also SPECEST_MTMFFT, SPECEST_TFR, SPECEST_HILBERT, SPECEST_MTMWELCH, SPECEST_NANFFT, SPECEST_MVAR, SPECEST_WLTCONVOL


% get the optional input arguments
keyvalcheck(varargin, 'optional', {'taper','pad','timeoi','timwin','freqoi','tapsmofrq'});
taper     = keyval('taper',       varargin); if isempty(taper),    taper   = 'dpss';     end
pad       = keyval('pad',         varargin); 
timeoi    = keyval('timeoi',      varargin); if isempty(timeoi),   timeoi  = 'all';      end
timwin    = keyval('timwin',      varargin); 
freqoi    = keyval('freqoi',      varargin); if isempty(freqoi),   freqoi  = 'all';      end
tapsmofrq = keyval('tapsmofrq',   varargin);

% throw errors for required input
if isempty(tapsmofrq) && strcmp(taper, 'dpss')
  error('you need to specify tapsmofrq when using dpss tapers')
end
if isempty(timwin)
  error('you need to specify timwin')
elseif (length(timwin) ~= length(freqoi) && ~strcmp(freqoi,'all'))
  error('timwin should be of equal length as freqoi')
end


% Set n's
[nchan,ndatsample] = size(dat);


% Determine fsample and set total time-length of data
fsample = 1/(time(2)-time(1));
dattime = ndatsample / fsample; % total time in seconds of input data

% Zero padding
if pad < dattime
  error('the padding that you specified is shorter than the data');
end
if isempty(pad) % if no padding is specified padding is equal to current data length
  pad = dattime;
end
postpad = zeros(1,round((pad - dattime) * fsample));
endnsample = pad * fsample;  % total number of samples of padded data
endtime    = pad;            % total time in seconds of padded data




% Set freqboi and freqoi
if isnumeric(freqoi) % if input is a vector
  freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1;
  freqoi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output
elseif strcmp(freqoi,'all') 
  freqboilim = round([0 fsample/2] ./ (fsample ./ endnsample)) + 1;
  freqboi    = freqboilim(1):1:freqboilim(2);
  freqoi     = (freqboi-1) ./ endtime;
end
% check for freqoi = 0 and remove it, there is no wavelet for freqoi = 0
if freqoi(1)==0
  freqoi(1)  = [];
  freqboi(1) = [];
  if length(timwin) == (length(freqoi) + 1)
    timwin(1) = [];
  end
end
nfreqboi = length(freqboi);
nfreqoi  = length(freqoi);



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


% set number of samples per time-window (timwin is in seconds)
timwinsample = round(timwin .* fsample);



% Compute tapers per frequency, multiply with wavelets and compute their fft
wltspctrm = cell(nfreqoi,1);
ntaper    = zeros(nfreqoi,1);
for ifreqoi = 1:nfreqoi
  
  switch taper
    case 'dpss'
      % create a sequence of DPSS tapers, ensure that the input arguments are double precision
      tap = double_dpss(timwinsample(ifreqoi), timwinsample(ifreqoi) .* (tapsmofrq(ifreqoi) ./ fsample))';
      % remove the last taper because the last slepian taper is always messy
      tap = tap(1:(end-1), :);

      % give error/warning about number of tapers
      if isempty(tap)
        error('datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',ndatsample/fsample,tapsmofrq(ifreqoi),fsample/fsample);
      elseif size(tap,1) == 1
        warning('using only one taper for specified smoothing')
      end
      
      
    case 'sine'
      tap = sine_taper(timwinsample(ifreqoi), timwinsample(ifreqoi) .* (tapsmofrq(ifreqoi) ./ fsample))';
      
    case 'alpha'
      tap = alpha_taper(timwinsample(ifreqoi), freqoi(ifreqoi)./ fsample)';
      tap = tap./norm(tap)';
      
    otherwise
      % create a single taper according to the window specification as a replacement for the DPSS (Slepian) sequence
      tap = window(taper, timwinsample(ifreqoi))';
      tap = tap ./ norm(tap,'fro'); % make it explicit that the frobenius norm is being used
  end
  
  % set number of tapers
  ntaper(ifreqoi) = size(tap,1);
  
  % Wavelet construction
  tappad   = ceil(endnsample ./ 2) - floor(timwinsample(ifreqoi) ./ 2);
  prezero  = zeros(1,tappad);
  postzero = zeros(1,round(endnsample) - ((tappad-1) + timwinsample(ifreqoi))-1);
  anglein  = (0:timwinsample(ifreqoi)-1)' .* ((2.*pi./fsample) .* freqoi(ifreqoi));
  wltspctrm{ifreqoi} = complex(zeros(size(tap,1),round(endnsample)));


  % the following code determines the phase-shift needed so that the centre of each wavelet has angle = 0. This code can probably be optimized greatly.
  % determine appropriate phase-shift so angle(wavelet) at center approximates 0 NOTE: this procedure becomes inaccurate when there are very few samples per cycle (i.e. 4-5)
  cyclefraction  = anglein / (2*pi); % transform angle to fraction of cycles
  if ((length(cyclefraction(cyclefraction<1))-1) < 5) % could be more robust
    warning('number of samples per wavelet cycle is less than 5')
  end
  fullcyclenum   = floor(max(cyclefraction)); % get the number of complete cycles in angle
  [dum fractind] = min(abs(cyclefraction - fullcyclenum)); % determine closest breakpoint in angle for which the last uncomplete cycle starts (closest so angle(wavelet) at centre gets closest to 0)
  if cyclefraction(fractind) < fullcyclenum % if index is from the last full cycle, shift it by 1. should be integrated in above line
    fractind = fractind + 1;
  end
  fractind = fractind + (length(prezero) - length(postzero)); % correct for unevend zero-padding (which shift the wavelet itself)
  fractind = (fractind + 1):length(anglein); % shift one sample upwards and fill indices (why again?)
  if length(fractind) > 1 % only continue if more than one sample can be split up
    % create new anglein with non-full cycly being split to both sides of the (resulting) wavelet
    nsplit     = length(fractind) / 2;
    anglestart = -(floor(nsplit):-1:1) .* ((2.*pi./fsample) .* freqoi(ifreqoi)); % using floor(nsplit) here, as the beginning of anglein is always at the exact start of sin/cos
    angleind   = 1:fractind(ceil(nsplit)); % using ceil(nsplit) here, as the end of anglein is nearly never at the end of a sin/cos, so an extra sample in case of non-integer-nsplit would do the most good here
    anglein    = [anglestart' ; anglein(angleind)];
  end
  
  for itap = 1:ntaper(ifreqoi)
    try % this try loop tries to fit the wavelet into wltspctrm, when its length is smaller than ndatsample, the rest is 'filled' with zeros because of above code
      % if a wavelet is longer than ndatsample, it doesn't fit and it is kept at zeros, which is translated to NaN's in the output
      % construct the complex wavelet
      coswav  = horzcat(prezero, tap(itap,:) .* cos(anglein)', postzero);
      sinwav  = horzcat(prezero, tap(itap,:) .* sin(anglein)', postzero);
      
      % consistency: cos must always be 1 at centre (necessary for angle(wavelet) at centre approximates 0), and sin must always be centered in upgoing flank (arbitrary), and
      centreind = round(length(coswav) / 2);
      % first the cos
      if coswav(centreind) < 0
        coswav = -coswav;
      end
      % now the sin
      if sinwav(centreind) > sinwav(centreind+1)
        sinwav = -sinwav;
      end
      wavelet = complex(coswav, sinwav);
      % debug plotting
      %figure; subplot(2,1,1);hold on;plot(real(wavelet));plot(imag(wavelet),'color','r'); tline = length(wavelet)/2;line([tline tline],[-0.2 0.2]); subplot(2,1,2);plot(angle(wavelet),'color','g');line([tline tline],[-pi pi])
      % store the fft of the complex wavelet
      wltspctrm{ifreqoi}(itap,:) = fft(wavelet,[],2);
    end
  end
end


% compute fft, major speed increases are possible here, depending on which matlab is being used whether or not it helps, which mainly focuses on orientation of the to be fft'd matrix
spectrum = complex(nan([sum(ntaper),nchan,nfreqoi,ntimeboi]));
datspectrum = fft([dat repmat(postpad,[nchan, 1])],[],2); 
for ifreqoi = 1:nfreqoi
  fprintf('processing frequency %d (%.2f Hz), %d tapers\n', ifreqoi,freqoi(ifreqoi),ntaper(ifreqoi));
  for itap = 1:ntaper(ifreqoi)
    for ichan = 1:nchan
      % compute indices that will be used to extracted the requested fft output    
      nsamplefreqoi    = timwin(ifreqoi) .* fsample;
      reqtimeboiind    = find((timeboi >=  (nsamplefreqoi ./ 2)) & (timeboi <    ndatsample - (nsamplefreqoi ./2)));
      reqtimeboi       = timeboi(reqtimeboiind);
      
      % compute datspectrum*wavelet, if there are reqtimeboi's that have data
      if ~isempty(reqtimeboi)
        dum = fftshift(ifft(datspectrum(ichan,:) .* wltspctrm{ifreqoi}(itap,:),[],2)); % fftshift is necessary because of post zero-padding, not necessary when pre-padding
        spectrum(itap,ichan,ifreqoi,reqtimeboiind) = dum(reqtimeboi);
      end
    end
  end
end




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for Matlab 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin)
tap = dpss(double(a), double(b), varargin{:});

















