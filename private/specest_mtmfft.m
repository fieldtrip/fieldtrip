function [spectrum,foi] = specest_mtmfft(dat, time, varargin) 

% SPECEST_MTMFFT computes a fast Fourier transform using many possible tapers
%
%
% Use as
%   [spectrum,foi] = specest_mtmfft(dat,time...)   
%
%   dat      = matrix of chan*sample 
%   time     = vector, containing time in seconds for each sample
%   spectrum = matrix of taper*chan*sample of fourier coefficients
%   foi      = vector of frequencies in spectrum
%
%
%
%
% Optional arguments should be specified in key-value pairs and can include:
%   taper      = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%   pad        = number, total length of data after zero padding (in seconds)
%   foi        = vector, containing frequencies of interest                                           
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
%

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'taper','pad','foi','tapsmofrq'});
taper     = keyval('taper',       varargin); if isempty(taper),    taper   = 'dpss';     end
pad       = keyval('pad',         varargin);
foi       = keyval('foi',         varargin); if isempty(foi),      foi     = 'max';      end  
tapsmofrq = keyval('tapsmofrq',   varargin); %%%% NOW CAN ONLY BE A NUMBER, IN MTMCONVOL IT CAN BE A VECTOR


% Set n's
[nchan,nsample] = size(dat);


% Determine fsample
fsample = nsample / (time(end) - time(1));


% Zero padding
if pad < (time(end) - time(1))
  error('the padding that you specified is shorter than the data');
end
if isempty(pad) % if no padding is specified padding is equal to current data length
  pad = (time(end)-time(1));
end
postpad = zeros(1,round(((pad - (time(end)-time(1))) * fsample) ./ 2)); % 'postpad', so naming concurs with mtmconvol



% Set fboi and foi 
if isnumeric(foi) % if input is a vector
  fboi    = round(foi ./ (fsample ./ (pad * fsample))) + 1;
  nfboi   = size(fboi,2);
  foi     = (fboi-1) ./ pad; % boi - 1 because 0 Hz is included in fourier output..... is this going correctly?
elseif strcmp(foi,'max') % if input was 'max'
  fboilim = round([0 fsample/2] ./ (fsample ./ (pad * fsample))) + 1;
  fboi    = fboilim(1):fboilim(2);
  nfboi   = size(fboi,2);
  foi     = (fboi-1) ./ pad;
end
if isempty(tapsmofrq) % default tapsmofrq
  tapsmofrq = 4;
end


% create tapers
switch taper
   
  case 'dpss'
    % create a sequence of DPSS tapers, ensure that the input arguments are double precision
    tap = double_dpss(nsample,nsample*(tapsmofrq./fsample))';
    % remove the last taper because the last slepian taper is always messy
    tap = tap(1:(end-1), :);
    
    % give error/warning about number of tapers
    if isempty(tap)
      error('datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',nsample/fsample,tapsmofrq,fsample/fsample);
    elseif size(tap,1) == 1
      warning('using only one taper for specified smoothing')
    end
        
  case 'sine'
    tap = sine_taper(nsample, nsample*(tapsmofrq./fsample))';
    
  case 'alpha'
    error('not yet implemented');
    
  otherwise
    % create the taper and ensure that it is normalized
    tap = window(taper, nsample)';
    tap = tap ./ norm(tap,'fro');
    
end % switch taper
ntap = size(tap,1);


% compute fft per channel, keeping tapers automatically (per channel is about 40% faster than all channels at the same time)
% compute fft, major speed increases are possible here, depending on which matlab is being used whether or not it helps, which mainly focuses on orientation of the to be fft'd matrix
spectrum = complex(zeros(ntap,nchan,nfboi),zeros(ntap,nchan,nfboi));
for itap = 1:ntap
  for ichan = 1:nchan
    dum = fft([dat(ichan,:) .* tap(itap,:) postpad],[],2); % would be much better if fft could take boi as input (muuuuuch less computation)
    spectrum(itap,ichan,:) = dum(fboi);
  end
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for Matlab 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin)
tap = dpss(double(a), double(b), varargin{:});


