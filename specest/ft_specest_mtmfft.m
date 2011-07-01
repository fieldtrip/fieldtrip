function [spectrum,ntaper,freqoi] = ft_specest_mtmfft(dat, time, varargin) 

% SPECEST_MTMFFT computes a fast Fourier transform using multitapering with
% the DPSS sequence or using a variety of single tapers
%
% Use as
%   [spectrum,freqoi] = specest_mtmfft(dat,time...)   
% where
%   dat      = matrix of chan*sample 
%   time     = vector, containing time in seconds for each sample
%   spectrum = matrix of taper*chan*freqoi of fourier coefficients
%   ntaper   = vector containing number of tapers per element of freqoi
%   freqoi   = vector of frequencies in spectrum
%
% Optional arguments should be specified in key-value pairs and can include:
%   taper      = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%   pad        = number, total length of data after zero padding (in seconds)
%   freqoi     = vector, containing frequencies of interest                                           
%   tapsmofrq  = the amount of spectral smoothing through multi-tapering. Note: 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box
%
% See also SPECEST_MTMCONVOL, SPECEST_CONVOL, SPECEST_HILBERT, SPECEST_WAVELET

% Copyright (C) 2010, Donders Institute for Brain, Cognition and Behaviour
%
% $Log$

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'taper','pad','freqoi','tapsmofrq','feedback','polyremoval'});
taper     = keyval('taper',       varargin); if isempty(taper),  error('You must specify a taper');    end
pad       = keyval('pad',         varargin);
freqoi    = keyval('freqoi',      varargin); if isempty(freqoi),   freqoi  = 'all';      end  
tapsmofrq = keyval('tapsmofrq',   varargin); 
fbopt     = keyval('feedback',    varargin);
polyorder = keyval('polyremoval', varargin); if isempty(polyorder), polyorder = 1; end

if isempty(fbopt),
  fbopt.i = 1;
  fbopt.n = 1;
end

% throw errors for required input
if isempty(tapsmofrq) && strcmp(taper, 'dpss')
  error('you need to specify tapsmofrq when using dpss tapers')
end

% Set n's
[nchan,ndatsample] = size(dat);

% Remove polynomial fit from the data -> default is demeaning
if polyorder >= 0
  dat = ft_preproc_polyremoval(dat, polyorder, 1, ndatsample);
end

% Determine fsample and set total time-length of data
fsample = 1/(time(2)-time(1));
dattime = ndatsample / fsample; % total time in seconds of input data

% Zero padding
if round(pad * fsample) < ndatsample
  error('the padding that you specified is shorter than the data');
end
if isempty(pad) % if no padding is specified padding is equal to current data length
  pad = dattime;
end
postpad    = zeros(1,ceil((pad - dattime) * fsample));
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;                   % total time in seconds of padded data

% Set freqboi and freqoi
if isnumeric(freqoi) % if input is a vector
  freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1;
  freqboi   = unique(freqboi);
  freqoi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output
elseif strcmp(freqoi,'all') % if input was 'all' 
  freqboilim = round([0 fsample/2] ./ (fsample ./ endnsample)) + 1;
  freqboi    = freqboilim(1):1:freqboilim(2);
  freqoi     = (freqboi-1) ./ endtime;
end
nfreqboi = length(freqboi);
nfreqoi  = length(freqoi);

% create tapers
switch taper
   
  case 'dpss'
    % create a sequence of DPSS tapers, ensure that the input arguments are double precision
    tap = double_dpss(ndatsample,ndatsample*(tapsmofrq./fsample))';
    % remove the last taper because the last slepian taper is always messy
    tap = tap(1:(end-1), :);
    
    % give error/warning about number of tapers
    if isempty(tap)
      error('datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',ndatsample/fsample,tapsmofrq,fsample/ndatsample);
    elseif size(tap,1) == 1
      warning_once('using only one taper for specified smoothing');
    end
        
  case 'sine'
    tap = sine_taper(ndatsample, ndatsample*(tapsmofrq./fsample))';
    tap = tap(1:(end-1), :); % remove the last taper
  
  case 'sine_old'
    % to provide compatibility with the tapers being scaled (which was default 
    % behavior prior to 29apr2011) yet this gave different magnitude of power 
    % when comparing with slepian multi tapers
    tap = sine_taper_scaled(ndatsample, ndatsample*(tapsmofrq./fsample))';
    tap = tap(1:(end-1), :); % remove the last taper
  
  case 'alpha'
    error('not yet implemented');
    
  case 'hanning'
    tap = hanning(ndatsample)'; 
    tap = tap./norm(tap, 'fro');

  otherwise
    % create the taper and ensure that it is normalized
    tap = window(taper, ndatsample)';
    tap = tap ./ norm(tap,'fro');
    
end % switch taper
ntaper = repmat(size(tap,1),nfreqoi,1);

str = sprintf('nfft: %d samples, datalength: %d samples, %d tapers',endnsample,ndatsample,ntaper(1));

[st, cws] = dbstack;
if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis')
  % specest_mtmfft has been called by ft_freqanalysis, meaning that ft_progress has been initialised
  ft_progress(fbopt.i./fbopt.n, ['processing trial %d/%d ',str,'\n'], fbopt.i, fbopt.n);
else
  fprintf([str, '\n']);
end

% determine phase-shift so that for all frequencies angle(t=0) = 0
timedelay = time(1); 
if timedelay ~= 0
  angletransform = complex(zeros(1,nfreqoi));
  for ifreqoi = 1:nfreqoi
    missedsamples = round(timedelay * fsample);
    % determine angle of freqoi if oscillation started at 0
    % the angle of wavelet(cos,sin) = 0 at the first point of a cycle, with sine being in upgoing flank, which is the same convention as in mtmconvol
    anglein = (missedsamples) .* ((2.*pi./fsample) .* freqoi(ifreqoi));
    coswav  = cos(anglein);
    sinwav  = sin(anglein);
    angletransform(ifreqoi) = angle(complex(coswav,sinwav));
  end
  angletransform = repmat(angletransform,[nchan,1]);
end

% compute fft, major speed increases are possible here, depending on which matlab is being used whether or not it helps, which mainly focuses on orientation of the to be fft'd matrix
spectrum = cell(ntaper(1),1);
for itap = 1:ntaper(1)
    dum = transpose(fft(transpose([dat .* repmat(tap(itap,:),[nchan, 1]) repmat(postpad,[nchan, 1])]))); % double explicit transpose to speedup fft
    dum = dum(:,freqboi);
    % phase-shift according to above angles
    if timedelay ~= 0
      dum = dum .* exp(-1i*angletransform);
    end
    dum = dum .* sqrt(2 ./ endnsample);
    spectrum{itap} = dum;
end
spectrum = reshape(vertcat(spectrum{:}),[nchan ntaper(1) nfreqboi]);% collecting in a cell-array and later reshaping provides significant speedups
spectrum = permute(spectrum, [2 1 3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for Matlab 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin)
tap = dpss(double(a), double(b), varargin{:});
