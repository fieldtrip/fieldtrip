function [spectrum,foi,toi] = specest_hilbert(dat, time, varargin)

% SPECEST_HILBERT performs a spectral estimation of data by repeatedly
% applying a bandpass filter and then doing a hilbert transform.
%
% Use as
%   [spectrum,foi,toi] = specest_hilbert(dat,time,...)  
%
%   dat      = matrix of chan*sample
%   time     = vector, containing time in seconds for each sample
%   spectrum = matrix of taper*chan*foi*toi of fourier coefficients
%   foi      = vector of frequencies in spectrum
%   toi      = vector of timebins in spectrum
%
%
%
%
% Optional arguments should be specified in key-value pairs and can include:
%   toi       = vector, containing time points of interest (in seconds)
%   foi       = vector, containing frequencies (in Hz)
%   width     = 
%   filttype  = 
%   filtorder = 
%   filtdir   = 
%
%
%
%
% See also SPECEST_MTMFFT, SPECEST_MTMCONVOL, SPECEST_TFR

% Copyright (C) 2010, Robert Oostenveld
%
% $Rev$

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'foi','foi','width','filttype','filtorder','filtdir'});
foi       = keyval('foi',       varargin);
toi       = keyval('toi',       varargin);   if isempty(toi),      toi     = 'max';      end
width     = keyval('width',     varargin);   if isempty(width),    width   = 1;          end
filttype  = keyval('filttype',  varargin);
filtorder = keyval('filtorder', varargin);
filtdir   = keyval('filtdir',   varargin);


% determine the data characteristics
[nchans, ntime] = size(dat);
fsample = 1/(time(2)-time(1));


% set a default sampling for the frequencies-of-interest
if isempty(foi),
  foi = linspace(2*width, (fsample/3), 50);
end
% check for foi = 0 and remove it
if any(foi==0)
  foi(foi==0) = [];
end
nfreq = length(foi);


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



% each frequency can have its own width
if numel(width)==1
  width = width*ones(size(foi));
end

% preallocate the result
spectrum = complex(zeros(nchans, nfreq, ntime));

for i=1:nfreq
  flt = preproc_bandpassfilter(dat, fsample, [foi(i)-width(i) foi(i)+width(i)], filtorder, filttype, filtdir);
  spectrum(:,i,:) = transpose(hilbert(transpose(flt)));
end

% get tboi out of spectrum
spectrum = spectrum(:,:,tboi);





