function [out,foi,toi] = specest_mtmconvol(dat, time, varargin)

% SPECEST_MTMCONVOL performs wavelet convolution in the time domain by multiplication in the frequency domain
%
%
% Use as
%   [out,foi,toi] = specest_mtmfft(dat,time,...)  %%% DECIDE: WHICH INPUT IS REQUIRED AND WHICH WILL BE DEFAULTED?
%
%   dat      = matrix of chan*sample
%   time     = vector, containing time in seconds for each sample
%   out      = matrix of taper*chan*sample of fourier coefficients
%   foi      = vector of frequencies in out
%   toi      = vector of timebins in out
%
%
%
%
% Optional arguments should be specified in key-value pairs and can include:
%   taper     = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%   pad       = number, indicating time-length of data to be padded out to          %%% IS IN SECONDS, MTMFFT HAS IT SAMPLES ATM
%   toi       = vector, containing time points of interest (in seconds, analysis window will be centered around these time points)
%   timwin    = vector, containing length of time windows (in seconds)
%   foi       = vector, containing frequencies (in Hz)
%   tapsmofrq = vector, the amount of spectral smoothing through multi-tapering. Note: 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box
%
%
%
%
%
%
%
%

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'dpss','pad','toi','timwin','foi','tapsmofrq'});
taper     = keyval('dpss',        varargin); if isempty(taper),    taper   = 'dpss';     end
pad       = keyval('pad',         varargin); 
toi       = keyval('toi',         varargin); if isempty(toi),      toi     = 'max';      end
timwin    = keyval('timwin',      varargin); % will be defaulted below
foi       = keyval('foi',         varargin); if isempty(foi),      foi     = 'max';      end
tapsmofrq = keyval('tapsmofrq',   varargin);

% Set n's
[nchan,nsample] = size(dat);

% Determine fsample
fsample = nsample / (time(end) - time(1));

% Zero padding (if pad is empty, no padding wil be performed in the end)
if pad < (time(end) - time(1))
  error('the padding that you specified is shorter than the data');
end
if isempty(pad) || (pad == (time(end) - time(1))) % if no padding is specified or if padding is equal to current data length
  pad = (time(end)-time(1));
else
  pad = round(((pad - (time(end) - time(1))) * fsample) / 2); 
end
prepad  = zeros(1,pad);
postpad = zeros(1,pad);



% Set fboi and foi and default tapsmofrq
% fboi = frequency bins of interest (index to be used in fft-output)
% foi  = frequencies of interest (Hz)
if isnumeric(foi) % if input is a vector
  fboi   = round(foi ./ pad) + 1;
  numboi = size(fboi,2);
  foi    = (fboi-1) ./ pad; % boi - 1 because 0 Hz is included in fourier output..... is this going correctly?
elseif strcmp(foi,'max') % if input was 'max'
  fboilim = round([0 fsample/2] ./ (pad)) + 1;
  fboi    = fboilim(1):1:fboilim(2);
  numboi = size(fboi,2);
  foi    = (fboi-1) ./ pad;
end
nfoi = length(foi);
if isempty(tapsmofrq) % default tapsmofrq
  tapsmofrq = ones(nfoi,1)*4;
end


% Set tboi and toi
% tboi = time bins of interest (index to be used in fft-output)
% toi  = time-points of interest (s)
if isnumeric(toi) % if input is a vector
  tboi = round(toi .* fsample) + 1;
  toi  = round(toi .* fsample) ./ fsample;
elseif strcmp(toi,'max') % if input was 'max'
  toi  = time;
  tboi = 1:length(time);
end


% time-windows %%% FIXME: IF 0 HZ IS REQUESTED, THIS BREAKS
if isempty(timwin)
  timwin = (1 ./ foi) * 3; % default is 3 cycles of a frequency
end
% set number of samples per time-window
timwinsmp = timwin * fsample;



% Compute tapers per frequency and compute their fft
wltspctrm = cell(nfoi,1);
ntaper    = zeros(nfoi,1);
for ifoi = 1:nfoi
  
  switch taper
    case 'dpss'
      % create a sequence of DPSS tapers, ensure that the input arguments are double precision
      tap = double_dpss(timwinsmp(ifoi), timwinsmp(ifoi) .* (tapsmofrq(ifoi) ./ fsample));
      % remove the last taper
      tap = tap(1:(end-1), :);
      
    case 'sine'
      tap = sine_taper(timwinsmp(ifoi), timwinsmp(ifoi) .* (tapsmofrq(ifoi) ./ fsample));
      
    case 'alpha'
      tap = alpha_taper(timwinsmp(ifoi), foi(ifoi)./ fsample);
      tap = tap./norm(tap);
      
    otherwise
      % create a single taper according to the window specification as a replacement for the DPSS (Slepian) sequence
      tap = window(cfg.taper, timwinsmp(ifoi));
      tap = tap ./ norm(tap);
  end
  
  ntaper(ifoi) = size(tap,2);
  % give error/warning about number of tapers
  if isempty(tap)
    error('datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',nsample/fsample,tapsmofrq,fsample/fsample);
  elseif size(tap,1) == 1
    warning('using only one taper for specified smoothing')
  end
  
  
  % Wavelet construction
  ins      = ceil(nsample ./ 2) - floor(timwinsmp(ifoi) ./ 2);
  prezero  = zeros(ins,1);
  postzero = zeros(nsample - ((ins-1) + timwinsmp(ifoi))-1,1);
  ind      = (0:timwinsmp(ifoi)-1)' .* ((2.*pi./fsample) .* foi(ifoi));
  wltspctrm{ifoi} = complex(zeros(size(tap,1),nsample));
  for itaper = 1:ntaper(ifoi)
    % construct the complex wavelet
    coswav  = vertcat(prezero, tap(:,itaper) .* cos(ind), postzero);
    sinwav  = vertcat(prezero, tap(:,itaper) .* sin(ind), postzero);
    wavelet = complex(coswav, sinwav);
    % store the fft of the complex wavelet
    wltspctrm{ifoi}(taplop,:) = fft(wavelet,[],1)';
  end
end






























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for Matlab 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin)
tap = dpss(double(a), double(b), varargin{:});

















