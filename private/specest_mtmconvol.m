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
%   taper      = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%   pad        = number, indicating time-length of data to be padded out to          %%% IS IN SECONDS, MTMFFT HAS IT SAMPLES ATM
%   toi       = vector, containing time points of interest (in seconds, analysis window will be centered around these time points)                                           
%   timwin     = vector, containing length of time windows (in seconds)
%   foi       = vector, containing frequencies (in Hz)                                              
%   tapsmofrq  = the amount of spectral smoothing through multi-tapering. Note: 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box
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
pad       = keyval('pad',         varargin); if isempty(pad),      pad     = [];         end 
toi       = keyval('toi',         varargin); if isempty(toi),      toi     = 'max';      end  
timwin    = keyval('timwin',      varargin); % will be defaulted below
foi       = keyval('foi',         varargin); if isempty(foi),      foi     = 'max';      end  
tapsmofrq = keyval('tapsmofrq',   varargin);

% set n's
[nchan,nsample] = size(dat);

% determine fsample
fsample = nsample / (time(end) - time(1));

% zero padding (if pad is empty, no padding wil be performed in the end)
if pad < (time(end) - time(1))
  error('the padding that you specified is shorter than the data');
end
pad = round(((pad - (time(end) - time(1))) * fsample) / 2);
prepad  = zeros(1,pad);
pstpad  = zeros(1,pad);


% set fboi and foi
if isnumeric(foi) % if input is a vector
fboi   = round(foi ./ (fsample ./ pad)) + 1;
numboi = size(fboi,2);
foi    = (fboi-1) ./ (pad / fsample); % boi - 1 because 0 Hz is included in fourier output..... is this going correctly?
elseif strcmp(foi,'max') % if input was 'max'
  fboilim = round([0 fsample/2] ./ (fsample ./ pad)) + 1;
  fboi    = fboilim(1):fboilim(2);
  numboi = size(fboi,2);
  foi    = (fboi-1) ./ (pad / fsample);
end

% set tboi and toi
if isnumeric(toi) % if input is a vector
elseif strcmp(toi,'max') % if input was 'max'
end






































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for Matlab 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin)
tap = dpss(double(a), double(b), varargin{:});

















