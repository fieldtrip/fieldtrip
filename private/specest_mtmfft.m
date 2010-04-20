function [spectrum,freqoi] = specest_mtmfft(dat, time, varargin) 

% SPECEST_MTMFFT computes a fast Fourier transform using many possible tapers
%
%
% Use as
%   [spectrum,freqoi] = specest_mtmfft(dat,time...)   
%
%   dat      = matrix of chan*sample 
%   time     = vector, containing time in seconds for each sample
%   spectrum = matrix of taper*chan*freqoi of fourier coefficients
%   freqoi   = vector of frequencies in spectrum
%
%
%
%
% Optional arguments should be specified in key-value pairs and can include:
%   taper      = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%   pad        = number, total length of data after zero padding (in seconds)
%   freqoi     = vector, containing frequencies of interest                                           
%   tapsmofrq  = the amount of spectral smoothing through multi-tapering. Note: 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box
%
%
%
%
%
%  TO DO:
%  - anchor phase to specific part of time-window (t=0 current suggestion), how to do it in the current format? need to do multiplication freq-domain (ala mtmconvol)?
%  - implement computation reduction by keeping tapers and such one way or another
%
%
% See also SPECEST_MTMCONVOL, SPECEST_HILBERT, SPECEST_TFR, SPECEST_NANFFT


% get the optional input arguments
keyvalcheck(varargin, 'optional', {'taper','pad','freqoi','tapsmofrq'});
taper     = keyval('taper',       varargin); if isempty(taper),  error('You must specify a taper');    end
pad       = keyval('pad',         varargin);
freqoi    = keyval('freqoi',      varargin); if isempty(freqoi),   freqoi  = 'max';      end  
tapsmofrq = keyval('tapsmofrq',   varargin); 

% throw errors for required input
if isempty(tapsmofrq) && strcmp(taper, 'dpss')
  error('you need to specify tapsmofrq when using dpss tapers')
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
postpad = zeros(1,ceil((pad - dattime) * fsample));
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;            % total time in seconds of padded data
postpad = [];


% Set freqboi and freqoi
if isnumeric(freqoi) % if input is a vector
  freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1;
  freqoi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output..... is this going correctly?
elseif strcmp(freqoi,'max') % if input was 'max' THIS IS IRRELEVANT, BECAUSE TIMWIN IS A REQUIRED INPUT NOW
  freqboilim = round([0 fsample/2] ./ (fsample ./ endnsample)) + 1;
  freqboi    = freqboilim(1):1:freqboilim(2);
  freqoi     = (freqboi-1) ./ endtime;
end
nfreqboi   = length(freqboi);
nfreqoi = length(freqoi);




% create tapers
switch taper
   
  case 'dpss'
    % create a sequence of DPSS tapers, ensure that the input arguments are double precision
    tap = double_dpss(ndatsample,ndatsample*(tapsmofrq./fsample))';
    % remove the last taper because the last slepian taper is always messy
    tap = tap(1:(end-1), :);
    
    % give error/warning about number of tapers
    if isempty(tap)
      error('datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',ndatsample/fsample,tapsmofrq,fsample/fsample);
    elseif size(tap,1) == 1
      warning('using only one taper for specified smoothing')
    end
        
  case 'sine'
    tap = sine_taper(ndatsample, ndatsample*(tapsmofrq./fsample))';
    
  case 'alpha'
    error('not yet implemented');
    
  otherwise
    % create the taper and ensure that it is normalized
    tap = window(taper, ndatsample)';
    tap = tap ./ norm(tap,'fro');
    
end % switch taper
ntap = size(tap,1);


% determine phase-shift so that for all frequencies angle(t=0) = 0
timedelay = abs(time(1)); % phase shift is equal for both negative and positive offsets
if timedelay ~= 0
  angletransform = complex(zeros(1,nfreqoi));
  for ifreqoi = 1:nfreqoi
    missedsamples = length(0:1/fsample:timedelay);
    % determine angle of freqoi if oscillation started at 0
    % the angle of wavelet(cos,sin) = 0 at the first point of a cycle, with sin being in upgoing flank, which is the same convention as in mtmconvol
    anglein = (missedsamples-1) .* ((2.*pi./fsample) .* freqoi(ifreqoi));
    coswav  = cos(anglein);
    sinwav  = sin(anglein);
    angletransform(ifreqoi) = angle(complex(coswav,sinwav));
  end
end

% compute fft per channel, keeping tapers automatically (per channel is about 40% faster than all channels at the same time)
% compute fft, major speed increases are possible here, depending on which matlab is being used whether or not it helps, which mainly focuses on orientation of the to be fft'd matrix
spectrum = complex(zeros(ntap,nchan,nfreqboi),zeros(ntap,nchan,nfreqboi));
for itap = 1:ntap
  for ichan = 1:nchan
    dum = fft([dat(ichan,:) .* tap(itap,:) postpad],[],2); % would be much better if fft could take boi as input (muuuuuch less computation)
    dum = dum(freqboi);
    % phase-shift according to above angles
    if timedelay ~= 0
      dum = dum .* (exp(-1i*(angle(dum) - angletransform)));
    end
    spectrum(itap,ichan,:) = dum;
  end
end
fprintf('nfft: %d samples, taper length: %d samples, %d tapers\n',endnsample,ndatsample,ntap);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for Matlab 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin)
tap = dpss(double(a), double(b), varargin{:});


