function [spectrum, freqoi, timeoi] = ft_specest_convol(dat, time, varargin)

% FT_SPECEST_CONVOL performs wavelet convolution in the time domain by
% convolution with Morlet's wavelets.
%
% Use as
%   [spectrum,freqoi,timeoi] = specest_convol(dat,time,...)
% where
%   dat      = matrix of chan*sample
%   time     = vector, containing time in seconds for each sample
%   spectrum = matrix of chan*freqoi*timeoi of fourier coefficients
%   freqoi   = vector of frequencies in spectrum
%   timeoi   = vector of timebins in spectrum
%
% Optional arguments should be specified in key-value pairs and can include
%   timeoi        = vector, containing time points of interest (in seconds, analysis window will be centered around these time points)
%   freqoi        = vector, containing frequencies (in Hz)
%   waveletwidth  = number, 'width' of wavelets expressed in cycles (default = 7)
%
% See also FT_FREQANALYSIS, FT_SPECEST_MTMFFT, FT_SPECEST_MTMCONVOL, FT_SPECEST_HILBERT, FT_SPECEST_NANFFT, FT_SPECEST_WAVELET

% Copyright (C) 2010, Donders Institute for Brain, Cognition and Behaviour
%
% $Log$

% get the optional input arguments
timeoi        = ft_getopt(varargin, 'timeoi', 'all');
freqoi        = ft_getopt(varargin, 'freqoi', 'all');
waveletwidth  = ft_getopt(varargin, 'waveletwidth', 7);
polyorder     = ft_getopt(varargin, 'polyorder', 1);

if isempty(fbopt),
  fbopt.i = 1;
  fbopt.n = 1;
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
endnsample = ndatsample;  % for consistency with mtmconvol and mtmfft
endtime    = dattime;     % for consistency with mtmconvol and mtmfft

% Set freqboi and freqoi
if isnumeric(freqoi) % if input is a vector
  freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1;
  freqboi   = unique(freqboi);
  freqoi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output
elseif strcmp(freqoi,'all')
  freqboilim = round([0 fsample/2] ./ (fsample ./ endnsample)) + 1;
  freqboi    = freqboilim(1):1:freqboilim(2);
  freqoi     = (freqboi-1) ./ endtime;
end
nfreqboi   = length(freqboi);
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

% compute wavelet family
wavfam = waveletfam(freqoi,fsample,waveletwidth);

% compute spectrum by convolving the wavelets with the data
spectrum = zeros(nchan, nfreqoi, nsample);
for ifreqoi = 1:nfreqoi
  str = sprintf('frequency %d (%.2f Hz)', ifreqoi,freqoi(ifreqoi));
  [st, cws] = dbstack;
  if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis') && verbose
    % specest_convol has been called by ft_freqanalysis, meaning that ft_progress has been initialised
    ft_progress(fbopt.i./fbopt.n, ['trial %d, ',str,'\n'], fbopt.i);
  elseif verbose
    fprintf([str, '\n']);
  end
  
  wavelet = wavfam{ifreqoi};
  for ichan = 1:nchan
    spectrum(ichan,ifreqoi,:) = conv(dat(ichan,:),  wavelet, 'same');
  end
  % pad the edges with nans to indicate that the wavelet was not fully immersed in the data THERE ARE NO NANS ADDED IN FREQANALYSIS_TFR?
  nanpad = ceil(length(wavfam{ifreqoi})/2);
  % the padding should not be longer than the actual data
  nanpad = min(nanpad, nsample);
  begnanpad = 1:nanpad;
  endnanpad = (nsample-nanpad+1):nsample;
  spectrum(:,ifreqoi,begnanpad) = nan;
  spectrum(:,ifreqoi,endnanpad) = nan;
end

% select the samples for the output
spectrum = spectrum(:,:,timeboi);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for waveletanalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = waveletfam(foi,fsample,waveletwidth)
dt = 1/fsample;
for k=1:length(foi)
  sf  = foi(k)/waveletwidth;
  st  = 1/(2*pi*sf);
  toi = -3.5*st:dt:3.5*st;
  A   = 1/sqrt(st*sqrt(pi));
  M{k}= A*exp(-toi.^2/(2*st^2)).*exp(i*2*pi*foi(k).*toi);
end

