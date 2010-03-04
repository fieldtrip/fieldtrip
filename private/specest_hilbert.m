function [spec, freq, time] = specest_hilbert(dat, time, varargin)

% SPECEST_HILBERT performs a spectral estimation of data by repeatedly
% applying a bandpass filter and then doing a hilbert transform.
%
% See also SPECEST_MTMFFT, SPECEST_MTMCONVOL, SPECEST_TFR

% Copyright (C) 2010, Robert Oostenveld
%
% $Rev$

% get the optional input arguments
foi       = keyval('foi',       varargin);
width     = keyval('width',     varargin); if isempty(width), width = 1; end
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
nfreq = length(foi);

% each frequency can have its own width
if numel(width)==1
  width = width*ones(size(foi));
end

% preallocate the result
spec = complex(zeros(nchans, nfreq, ntime));
% these are also returned
freq = foi;
time = time;

for i=1:nfreq
  flt = preproc_bandpassfilter(dat, fsample, [foi(i)-width(i) foi(i)+width(i)], filtorder, filttype, filtdir);
  spec(:,i,:) = transpose(hilbert(transpose(flt)));
end
