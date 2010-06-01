function [spectrum,freqoi,timeoi] = specest_hilbert(dat, time, varargin)

% SPECEST_HILBERT performs a spectral estimation of data by repeatedly
% applying a bandpass filter and then doing a hilbert transform.
%
% Use as
%   [spectrum,freqoi,timeoi] = specest_hilbert(dat,time,...)  
%
%   dat      = matrix of chan*sample
%   time     = vector, containing time in seconds for each sample
%   spectrum = matrix of taper*chan*freqoi*timeoi of fourier coefficients
%   freqoi   = vector of frequencies in spectrum
%   timeoi   = vector of timebins in spectrum
%
%
%
%
% Optional arguments should be specified in key-value pairs and can include:
%   timeoi    = vector, containing time points of interest (in seconds)
%   freqoi    = vector, containing frequencies (in Hz)
%   width     = 
%   filttype  = 
%   filtorder = 
%   filtdir   = 
%
%
%
%
% See also SPECEST_MTMFFT, SPECEST_TFR, SPECEST_MTMCONVOL, SPECEST_MTMWELCH, SPECEST_NANFFT, SPECEST_MVAR, SPECEST_WLTCONVOL

% Copyright (C) 2010, Robert Oostenveld
%
% $Rev$

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'freqoi','timeoi','width','filttype','filtorder','filtdir'});
freqoi    = keyval('freqoi',    varargin);
timeoi    = keyval('timeoi',    varargin);   if isempty(timeoi),   timeoi  = 'all';      end
width     = keyval('width',     varargin);   if isempty(width),    width   = 1;          end
filttype  = keyval('filttype',  varargin);
filtorder = keyval('filtorder', varargin);
filtdir   = keyval('filtdir',   varargin);


% Set n's
[nchan,ndatsample] = size(dat);


% Determine fsample and set total time-length of data
fsample = 1/(time(2)-time(1));


% set a default sampling for the frequencies-of-interest
if isempty(freqoi),
  freqoi = linspace(2*width, (fsample/3), 50);
end
% check for freqoi = 0 and remove it
if any(freqoi==0)
  freqoi(freqoi==0) = [];
end
nfreq = length(freqoi);


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


% each frequency can have its own width
if numel(width)==1
  width = width*ones(size(freqoi));
end

% preallocate the result
spectrum = complex(zeros(nchans, nfreq, ntime));

for i=1:nfreq
  flt = preproc_bandpassfilter(dat, fsample, [freqoi(i)-width(i) freqoi(i)+width(i)], filtorder, filttype, filtdir);
  spectrum(:,i,:) = transpose(hilbert(transpose(flt)));
end

% get timeboi out of spectrum
spectrum = spectrum(:,:,timeboi);





