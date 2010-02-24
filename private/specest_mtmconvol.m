function [out,foi,toi] = specest_mtmconvol(dat, fsample, varargin) 

% SPECEST_MTMCONVOL performs wavelet convolution in the time domain by multiplication in the frequency domain
%
%
% Use as
%   [out,foi,toi] = specest_mtmfft(dat,fsample,...)  %%% DECIDE: WHICH INPUT IS REQUIRED AND WHICH WILL BE DEFAULTED?
%
%   dat      = matrix of chan*sample 
%   fsample  = number indicating sampling rate
%   out      = matrix of taper*chan*sample of fourier coefficients
%   foi      = vector of frequencies in out
%   toi      = vector of timebins in out
%
%
%
%
% Optional arguments should be specified in key-value pairs and can include:
%   taper      = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%   pad        = number, total number of samples after zero padding               %%% COULD ALSO BE GIVEN IN SECONDS....
%   time       = vector, containing time points of interest (in seconds, analysis window will be centered around these time points)                                           
%   timwin     = vector, containing length of time windows (in seconds)
%   freq       = vector, containing frequencies (in Hz)                                              
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
keyvalcheck(varargin, 'optional', {'dpss','pad','time','timwin','freq','tapsmofrq'});
taper     = keyval('dpss',        varargin); if isempty(taper),    taper   = 'dpss';     end
pad       = keyval('pad',         varargin); if isempty(pad),      pad     = 0;          end
time      = keyval('time',        varargin); if isempty(time),     time    = 'max';      end  
timwin    = keyval('timwin',      varargin); % will be defaulted below
freq      = keyval('freq',        varargin); if isempty(freq),     freq    = 'max';      end  
tapsmofrq = keyval('tapsmofrq',   varargin);














































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for Matlab 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin)
tap = dpss(double(a), double(b), varargin{:});

















