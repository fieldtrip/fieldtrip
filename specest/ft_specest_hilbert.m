function [spectrum,freqoi,timeoi] = ft_specest_hilbert(dat, time, varargin)

% FT_SPECEST_HILBERT performs a spectral estimation of data by repeatedly
% applying a bandpass filter and then doing a hilbert transform.
%
% Use as
%   [spectrum,freqoi,timeoi] = specest_hilbert(dat,time,...)
% where
%   dat      = matrix of chan*sample
%   time     = vector, containing time in seconds for each sample
%   spectrum = matrix of chan*freqoi*timeoi of fourier coefficients
%   freqoi   = vector of frequencies in spectrum
%   timeoi   = vector of timebins in spectrum
%
% Optional arguments should be specified in key-value pairs and can include:
%   timeoi    = vector, containing time points of interest (in seconds)
%   freqoi    = vector, containing frequencies (in Hz)
%   pad       = number, indicating time-length of data to be padded out to in seconds
%   width     =
%   filttype  =
%   filtorder =
%   filtdir   =
%
% See also FT_FREQANALYSIS, FT_SPECEST_MTMFFT, FT_SPECEST_CONVOL, FT_SPECEST_MTMCONVOL, FT_SPECEST_WAVELET

% Copyright (C) 2010, Robert Oostenveld
%
% $Log: 3162 $

% get the optional input arguments
freqoi    = ft_getopt(varargin, 'freqoi');
timeoi    = ft_getopt(varargin, 'timeoi', 'all');
width     = ft_getopt(varargin, 'width', 1);
filttype  = ft_getopt(varargin, 'filttype');    if isempty(filttype),  error('you need to specify filter type'),         end
filtorder = ft_getopt(varargin, 'filtorder');   if isempty(filtorder), error('you need to specify filter order'),        end
filtdir   = ft_getopt(varargin, 'filtdir');     if isempty(filtdir),   error('you need to specify filter direction'),    end
pad       = ft_getopt(varargin, 'pad');
polyorder = ft_getopt(varargin, 'polyorder', 1);
verbose   = ft_getopt(varargin, 'verbose', true);

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
prepad  = zeros(1,floor(((pad - dattime) * fsample)./2));
postpad = zeros(1,ceil(((pad - dattime) * fsample)./2));

% set a default sampling for the frequencies-of-interest
if isempty(freqoi),
  freqoi = linspace(2*width, (fsample/3), 50);
end
% check for freqoi = 0 and remove it
if any(freqoi==0)
  freqoi(freqoi==0) = [];
end
nfreqoi = length(freqoi);

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

% expand width to array if constant width
if numel(width) == 1
  width = ones(1,nfreqoi) * width;
end

% create filter frequencies and check validity
filtfreq = [];
invalidind = [];
for ifreqoi = 1:nfreqoi
  tmpfreq = [freqoi(ifreqoi)+width(ifreqoi) freqoi(ifreqoi)-width(ifreqoi)];
  if all((sign(tmpfreq) == 1))
    filtfreq(end+1,:) = tmpfreq;
  else
      invalidind = [invalidind ifreqoi];
      warning_once(sprintf('frequency %.2f Hz cannot be estimated with resolution %.2f Hz', freqoi(ifreqoi), width(ifreqoi)));
  end
end

freqoi(invalidind) = [];
nfreqoi = size(filtfreq,1);

% preallocate the result and perform the transform
spectrum = complex(nan(nchan, nfreqoi, ntimeboi), nan(nchan, nfreqoi, ntimeboi));
for ifreqoi = 1:nfreqoi
  if verbose
    fprintf('processing frequency %d (%.2f Hz)\n', ifreqoi,freqoi(ifreqoi));
  end
  % filter
  flt = ft_preproc_bandpassfilter(dat, fsample, filtfreq(ifreqoi,:), filtorder, filttype, filtdir);
  
  % transform and insert
  dum = transpose(hilbert(transpose([repmat(prepad,[nchan, 1]) flt repmat(postpad,[nchan, 1])])));
  spectrum(:,ifreqoi,:) = dum(:,timeboi);
end
