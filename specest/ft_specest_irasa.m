function [spectrum,ntaper,freqoi] = ft_specest_irasa(dat, time, varargin)

% FT_SPECEST_IRASA estimates the powerspectral arrythmic component of the 
% time-domain using Irregular-Resampling Auto-Spectral Analysis (IRASA)
%
% Use as
%   [spectrum,ntaper,freqoi] = ft_specest_irasa(dat,time...)
% where
%   dat        = matrix of chan*sample
%   time       = vector, containing time in seconds for each sample
%   spectrum   = matrix of taper*chan*freqoi of fourier coefficients
%   ntaper     = vector containing number of tapers per element of freqoi
%   freqoi     = vector of frequencies in spectrum
%
% Optional arguments should be specified in key-value pairs and can include
%   taper      = 'dpss', 'hanning' or many others, see WINDOW (default = 'hanning')
%   pad        = number, total length of data after zero padding (in seconds)
%   padtype    = string, indicating type of padding to be used (see ft_preproc_padding, default: zero)
%   freqoi     = vector, containing frequencies of interest
%   tapsmofrq  = the amount of spectral smoothing through multi-tapering. Note: 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box
%   dimord     = 'tap_chan_freq' (default) or 'chan_time_freqtap' for memory efficiency (only used when variable number slepian tapers)
%   polyorder  = number, the order of the polynomial to fitted to and removed from the data prior to the fourier transform (default = 0 -> remove DC-component)
%   taperopt   = additional taper options to be used in the WINDOW function, see WINDOW
%   verbose    = output progress to console (0 or 1, default 1)
%
% This implements: Wen H, Liu Z. Separating fractal and oscillatory components in the power spectrum of neurophysiological signal. Brain Topogr. 2016 Jan;29(1):13-26.
%   For application, see Stolk et al., Electrocorticographic dissociation of 
%   alpha and beta rhythmic activity in the human sensorimotor system. It
%   is recommended the user first sub-segments the data using ft_redefinetrial 
%   and specifies cfg.pad = 'nextpow2' when calling ft_frequencyanalysis in 
%   order to implement steps A and B of the original algorithm in Wen & liu.
%
% See also FT_FREQANALYSIS, FT_SPECEST_MTMFFT, FT_SPECEST_MTMCONVOL, FT_SPECEST_TFR, FT_SPECEST_HILBERT, FT_SPECEST_WAVELET

% Copyright (C) 2019, Arjen Stolk
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% these are for speeding up computation of tapers on subsequent calls
persistent previous_argin previous_tap

% get the optional input arguments
taper     = ft_getopt(varargin, 'taper'); if isempty(taper), ft_error('You must specify a taper'); end
pad       = ft_getopt(varargin, 'pad');
padtype   = ft_getopt(varargin, 'padtype', 'zero');
freqoi    = ft_getopt(varargin, 'freqoi', 'all');
tapsmofrq = ft_getopt(varargin, 'tapsmofrq');
dimord    = ft_getopt(varargin, 'dimord', 'tap_chan_freq');
fbopt     = ft_getopt(varargin, 'feedback');
verbose   = ft_getopt(varargin, 'verbose', true);
polyorder = ft_getopt(varargin, 'polyorder', 0);
tapopt    = ft_getopt(varargin, 'taperopt');
hset      = ft_getopt(varargin, 'hset', 1.1:0.05:1.9); % IRASA resampling factors

if isempty(fbopt)
  fbopt.i = 1;
  fbopt.n = 1;
end

% throw errors for required input
if isempty(tapsmofrq) && (strcmp(taper, 'dpss') || strcmp(taper, 'sine'))
  ft_error('you need to specify tapsmofrq when using dpss or sine tapers')
end

% this does not work on integer data
dat = cast(dat, 'double');

% Set n's
[nchan,ndatsample] = size(dat);
nhset = length(hset);

% This does not work on integer data
if ~isa(dat, 'double') && ~isa(dat, 'single')
  dat = cast(dat, 'double');
end

% Remove polynomial fit from the data -> default is demeaning
if polyorder >= 0
  dat = ft_preproc_polyremoval(dat, polyorder, 1, ndatsample);
end

% Determine fsample and set total time-length of data
fsample = 1./mean(diff(time));
dattime = ndatsample / fsample; % total time in seconds of input data

% Zero padding
if round(pad * fsample) < ndatsample
  ft_error('the padding that you specified is shorter than the data');
end
if isempty(pad) % if no padding is specified padding is equal to current data length
  pad = dattime;
end
postpad    = ceil((pad - dattime) * fsample);
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;                   % total time in seconds of padded data

% Set freqboi and freqoi
freqoiinput = freqoi;
if isnumeric(freqoi) % if input is a vector
  freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1; % is equivalent to: round(freqoi .* endtime) + 1;
  freqboi   = unique(freqboi);
  freqoi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output
elseif strcmp(freqoi,'all') % if input was 'all'
  freqboilim = round([0 fsample/2] ./ (fsample ./ endnsample)) + 1;
  freqboi    = freqboilim(1):1:freqboilim(2);
  freqoi     = (freqboi-1) ./ endtime;
end
nfreqboi = length(freqboi);
nfreqoi  = length(freqoi);
if (strcmp(taper, 'dpss') || strcmp(taper, 'sine')) && numel(tapsmofrq)~=1 && (numel(tapsmofrq)~=nfreqoi)
  ft_error('tapsmofrq needs to contain a smoothing parameter for every frequency when requesting variable number of slepian tapers')
end

% throw a warning if input freqoi is different from output freqoi
if isnumeric(freqoiinput)
  % check whether padding is appropriate for the requested frequency resolution
  rayl = 1/endtime;
  if any(rem(freqoiinput,rayl)) % not always the case when they mismatch
    ft_warning('padding not sufficient for requested frequency resolution, for more information please see the FAQs on www.ru.nl/neuroimaging/fieldtrip');
  end
  if numel(freqoiinput) ~= numel(freqoi) % freqoi will not contain double frequency bins when requested
    ft_warning('output frequencies are different from input frequencies, multiples of the same bin were requested but not given');
  else
    if any(abs(freqoiinput-freqoi) >= eps*1e6)
      ft_warning('output frequencies are different from input frequencies');
    end
  end
end

% determine whether tapers need to be recomputed
current_argin = {time, postpad, taper, tapsmofrq, freqoi, tapopt, dimord}; % reasoning: if time and postpad are equal, it's the same length trial, if the rest is equal then the requested output is equal
if isequal(current_argin, previous_argin)
  % don't recompute tapers
  tap = previous_tap;
else
  % (re)compute tapers, 1:mid are upsample tapers, mid+1:end are downsample tapers
  for ih = 1:nhset
    [n, d] = rat(hset(ih)); % n > d
    udat = resample(dat', n, d)'; % upsample
    utap = hanning(size(udat,2))'; % Hanning taper
    tap{ih,1} = utap./norm(utap, 'fro');
    ddat = resample(dat', d, n)'; % downsample
    dtap = hanning(size(ddat,2))'; % Hanning taper
    tap{ih+nhset,1} = dtap./norm(dtap, 'fro');
  end
end

% set ntaper
if ~((strcmp(taper,'dpss') || strcmp(taper,'sine')) && numel(tapsmofrq)>1) % variable number of slepian tapers not requested
  ntaper = repmat(size(tap,2),nfreqoi,1); % pretend there's only one taper
else % variable number of slepian tapers requested
  ntaper = cellfun(@size,tap,repmat({1},[1 nfreqoi]));
end

% compute fft
str = sprintf('nfft: %d samples, datalength: %d samples, %d tapers',endnsample,ndatsample,ntaper(1));
[st, cws] = dbstack;
if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis')
  % specest_mtmfft has been called by ft_freqanalysis, meaning that ft_progress has been initialised
  ft_progress(fbopt.i./fbopt.n, ['processing trial %d/%d ',str,'\n'], fbopt.i, fbopt.n);
elseif verbose
  fprintf([str, '\n']);
end
spectrum = cell(ntaper(1),1);
for itap = 1:ntaper(1)
  %%%% IRASA STARTS %%%%
  pow = zeros(nchan,nfreqboi,nhset);
  for ih = 1:nhset % loop across resampling factors
    % resample
    [n, d] = rat(hset(ih)); % n > d
    udat = resample(dat', n, d)'; % upsample
    ddat = resample(dat', d, n)'; % downsample
    
    % fft of upsampled data
    ucom = fft(ft_preproc_padding(bsxfun(@times,udat,tap{ih,1}), padtype, 0, postpad), endnsample, 2); % fft
    ucom = ucom(:,freqboi);
    ucom = ucom .* sqrt(2 ./ size(udat,2));
    
    % fft of downsampled data
    dcom = fft(ft_preproc_padding(bsxfun(@times,ddat,tap{ih+nhset,1}), padtype, 0, postpad), endnsample, 2); % fft
    dcom = dcom(:,freqboi);
    dcom = dcom .* sqrt(2 ./ size(ddat,2));
    
    % geometric mean for this resampling factor
    upow = abs(ucom).^2;
    dpow = abs(dcom).^2;
    pow(:,:,ih) = sqrt(upow.*dpow); 
  end
  
  % median across resampling factors
  spectrum{itap} = median(pow,3);
  %%%% IRASA ENDS %%%%
end

spectrum = reshape(vertcat(spectrum{:}),[nchan ntaper(1) nfreqboi]); % collecting in a cell-array and later reshaping provides significant speedups
spectrum = permute(spectrum, [2 1 3]);

% remember the current input arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
previous_argin = current_argin;
previous_tap   = tap;
