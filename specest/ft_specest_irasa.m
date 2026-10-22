function [spectrum, ntaper, freqoi, hset] = ft_specest_irasa(dat, time, varargin)

% FT_SPECEST_IRASA separates the fractal components from the orginal power spectrum
% using Irregular-Resampling Auto-Spectral Analysis (IRASA)
%
% Use as
%   [spectrum, ntaper, freqoi] = ft_specest_irasa(dat, time, ...)
% where the input arguments are
%   dat        = matrix of chan*sample
%   time       = vector, containing time in seconds for each sample
% and the output arguments are
%   spectrum   = matrix of taper*chan*freqoi of fourier coefficients
%   ntaper     = vector containing number of tapers per element of freqoi
%   freqoi     = vector of frequencies in spectrum
%
% Optional arguments should be specified in key-value pairs and can include
%   freqoi     = vector, containing frequencies of interest
%   output     = string, indicating type of output ('fractal' or 'original', default 'fractal')
%   pad        = number, total length of data after zero padding (in seconds)
%   padtype    = string, indicating type of padding to be used, can be 'zero', 'mean', 'localmean', 'edge', or 'mirror' (default = 'zero')
%   polyorder  = number, the order of the polynomial to fitted to and removed from the data prior to the Fourier transform (default = 0, which removes the DC-component)
%   verbose    = boolean, output progress to console (0 or 1, default 1)
%
% This implements: Wen.H. & Liu.Z.(2016), Separating fractal and oscillatory components in the power
% spectrum of neurophysiological signal. Brain Topogr. 29(1):13-26. The source code accompanying the
% original paper is avaible from https://purr.purdue.edu/publications/1987/1
%
% For more information about the difference between the current and previous version and how to use this
% function, please see https://www.fieldtriptoolbox.org/example/irasa/
%
% See also FT_FREQANALYSIS, FT_SPECEST_MTMFFT, FT_SPECEST_MTMCONVOL, FT_SPECEST_TFR, FT_SPECEST_HILBERT, FT_SPECEST_WAVELET

% undocumented options
%   taper     = string, type of taper to use, default 'hanning'
%   tapopt    = additional taper options, if the corresponding taper
%                 function has additional options 
%   tapsmofrq = scalar, smoothing parameter (when taper = 'dpss')
%   weightopt = string, adaptive weighting method (when taper = 'dpss'), default 'mean'
%   ntaper    = scalar, number of tapers to use (when taper = 'dpss'), influences
%                 bias/variance trade off
%   mfunc     = string, specifying how to deal with the individual
%                 resampled fft's, default 'median'. (median_nomean,
%                 none_mean, none_individual, median_dpss, trimmean)

% Copyright (C) 2026, Jan-Mathijs Schoffelen
% Copyright (C) 2019-2020, Rui Liu, Arjen Stolk
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
persistent previous_argin previous_tap previous_lambda

% get the optional input arguments
freqoi    = ft_getopt(varargin, 'freqoi', 'all');
output    = ft_getopt(varargin, 'output', 'fractal');
pad       = ft_getopt(varargin, 'pad');
padtype   = ft_getopt(varargin, 'padtype', 'zero');
polyorder = ft_getopt(varargin, 'polyorder', 1);
verbose   = ft_getopt(varargin, 'verbose', true);
fbopt     = ft_getopt(varargin, 'feedback');
hset      = ft_getopt(varargin, 'hset', []); % array of stretch/compression fractions
nwindow   = ft_getopt(varargin, 'nwindow', 10);
windowlength = ft_getopt(varargin, 'windowlength', 'auto');
taper     = ft_getopt(varargin, 'taper', 'hanning');
mfunc     = ft_getopt(varargin, 'mfunc', 'median');
tapopt    = ft_getopt(varargin, 'taperopt');
tapsmofrq = ft_getopt(varargin, 'tapsmofrq');
weightopt = ft_getopt(varargin, 'weightopt', 'mean');
ntaper    = ft_getopt(varargin, 'ntaper', []);

if strcmp(taper, 'dpss')
  switch weightopt
    case 'mean'
      adaptflag = 0;
    case 'eig'
      adaptflag = 1;
    case 'adapt'
      adaptflag = 2;
    case 'qiinv'
      %adaptflag = 3;
      ft_error('weightopt ''qiinv'' specified for dpss-tapering is currently not supported');
    otherwise
      ft_error('unknown weightopt specified for dpss-tapering');
  end
elseif ~strcmp(weightopt, 'mean')
  ft_warning('no dpss tapers are specified, so the value %s of weightopt does not have an effect', weightopt);
end

% the original implementation uses a windowed (pwelch like) estimation technique,
% where - per epoch - the data are chunked into 10 sub-windows, with a length of
% 'nextpow2' times of (0.9*nsmp) samples. 

verbose = istrue(verbose); % if the calling function has 'yes'/'no'/etc

% check output option
if ~strcmp(output, {'fractal','original'})
  ft_error('The current version ft_specest_irasa outputs ''fractal'' or ''original'' power only. For more information about the update, see https://www.fieldtriptoolbox.org/example/irasa/');
end

% this does not work on data that is not single or double precision
if ~isa(dat, 'double')
  ft_info('Casting data to double precision')
  dat = cast(dat, 'double');
end

% size of the data, and sampling frequency
[nchan, nsmp] = size(dat);
fsample       = 1./mean(diff(time));

% parameters needed for the pwelch-type of sliding window
if isequal(windowlength, 'auto')
  windownsample = 2^floor(log2(nsmp*0.9));
elseif isequal(windowlength, 'all')
  windownsample = nsmp;
  nwindow = 1;
else
  windownsample = round(fsample*windowlength);
end
subset_dist = floor((nsmp - windownsample)/(nwindow-1));
if nwindow==1
  subset_dist = 0;
end

% resampling ratio can be passed as option, default is defined above
if isempty(hset)
  hset = (1.1:0.05:1.9);
end
nhset = length(hset);

% remove polynomial fit from the data -> default is demeaning
if polyorder >= 0
  dat = ft_preproc_polyremoval(dat, polyorder, 1, nsmp);
end

% check whether zero padding is needed
if isempty(pad)
  % if no padding is specified this is set equal to the current data length
  pad = 2*windownsample/fsample;
end
if round(pad * fsample) < windownsample
  ft_error('the padding that you specified is shorter than the data');
end
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;                   % total time in seconds of padded data

% set freqboi and freqoi
freqoiinput = freqoi;
if isnumeric(freqoi) % if input is a vector
  freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1; % is equivalent to: round(freqoi .* endtime) + 1;
  freqboi   = unique(freqboi);
  freqoi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output
elseif strcmp(freqoi, 'all') % if input was 'all'
  freqboilim = round([0 fsample/2] ./ (fsample ./ endnsample)) + 1;
  freqboi    = freqboilim(1):1:freqboilim(2);
  freqoi     = (freqboi-1) ./ endtime;
end
nfreqboi = length(freqboi);
nfreqoi  = length(freqoi);

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


% determine whether tapers need to be computed
current_argin = {output, time, endnsample, freqoi, taper, hset, tapsmofrq, nwindow, ntaper}; % reasoning: if time and endnsample are equal, it's the same length trial, if the rest is equal then the requested output is equal
if isequal(current_argin, previous_argin)
  % don't recompute tapers
  tap = previous_tap;
  lambda = previous_lambda;
else
  if strcmp(output, 'fractal')
    % (re)compute tapers, 1:mid are upsample tapers, mid+1:end are downsample tapers
    tap = cell(nhset*2,1);
    for ih = 1:nhset
      [n, d] = rat(hset(ih)); % n > d

      nsmp_up = size(resample(zeros(windownsample,1), n, d),1);
      switch taper
        case 'hanning'
          tmp = hanning(nsmp_up)';
          tmp = tmp./norm(tmp, 'fro');
          tmplambda = 1; % not needed
        case 'dpss'
          nw = nsmp*(tapsmofrq./fsample);
          [tmp, tmplambda] = dpss(nsmp_up, nw);
          if isempty(ntaper)
            ntaper = size(tmp,2)-1;
          else
            ntaper = min(ntaper, size(tmp,2)-1);
          end
          tmp = tmp(:,1:ntaper)';
          tmplambda = tmplambda(1:ntaper);

        case {'sine' 'sine_old' 'alpha'}
          ft_error('taper = %s is not implemented in ft_specest_irasa', taper);
        otherwise
          % create the taper and ensure that it is normalized
          if isempty(tapopt) % some windowing functions don't support nargin>1, and window.m doesn't check it
            tmp = window(taper, nsmp_up)';
            tmp = tmp./norm(tmp, 'fro');
          else
            tmp = window(taper, nsmp_up, tapopt)';
            tmp = tmp./norm(tmp, 'fro');
          end
          tmplambda = 1;
      end
      tap{ih,1}    = tmp;% for upsampled subsets
      lambda{ih,1} = tmplambda;

      nsmp_down = size(resample(zeros(windownsample,1), d, n),1);
      switch taper  
        case 'hanning'
          tmp = hanning(nsmp_down)';
          tmp = tmp./norm(tmp, 'fro');
        case 'dpss'
          nw  = windownsample*(tapsmofrq./fsample);
          [tmp, tmplambda] = dpss(nsmp_down, nw);
          tmp = tmp(:,1:ntaper)';
          tmplambda = tmplambda(1:ntaper);

        case {'sine' 'sine_old' 'alpha'}
          ft_error('taper = %s is not implemented in ft_specest_irasa', taper);
        otherwise
          % create the taper and ensure that it is normalized
          if isempty(tapopt) % some windowing functions don't support nargin>1, and window.m doesn't check it
            tmp = window(taper, nsmp_down)';
            tmp = tmp./norm(tmp, 'fro');
          else
            tmp = window(taper, nsmp_down, tapopt)';
            tmp = tmp./norm(tmp, 'fro');
          end
      end
      tap{ih+nhset,1}    = tmp;% for downsampled subsets
      lambda{ih+nhset,1} = tmplambda;
    end
  elseif strcmp(output, 'original')
    switch taper
      case 'hanning'
        tap = hanning(windownsample)';
        tap = tap./norm(tap, 'fro');
        lambda = 1; % not needed

      case 'dpss'
        nw  = windownsample*(tapsmofrq./fsample);
        [tmp, lambda] = dpss(windownsample, nw);
        tap = tmp(:,1:ntaper)';
        lambda = lambda(1:ntaper);
      case {'sine' 'sine_old' 'alpha'}
        ft_error('taper = %s is not implemented in ft_specest_irasa', taper);
      otherwise
        % create the taper and ensure that it is normalized
        if isempty(tapopt) % some windowing functions don't support nargin>1, and window.m doesn't check it
          tap = window(taper, windownsample)';
          tap = tap./norm(tap, 'fro');
        else
          tap = window(taper, windownsample, tapopt)';
          tap = tap./norm(tap, 'fro');
        end
        lambda = 1;
    end
    tap = {tap};
    lambda = {lambda};
  end
end

% set ntaper
ntaper = size(tap{1},1);

% feedback of computing progress
if isempty(fbopt)
  fbopt.i = 1;
  fbopt.n = 1;
end
str = sprintf('nfft: %d samples, datalength: %d samples, %d tapers',endnsample,nsmp,ntaper);
st  = dbstack;
if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis')
  % ft_specest_irasa has been called by ft_freqanalysis, meaning that ft_progress has been initialised
  ft_progress(fbopt.i./fbopt.n, ['processing trial %d/%d ',str,'\n'], fbopt.i, fbopt.n);
elseif verbose
  fprintf([str, '\n']);
end

%%%%%%%% IRASA %%%%%%%%%%
if strcmp(output,'fractal')
  if isequal(mfunc, 'none_individual')
    pow  = zeros(nchan,nfreqboi,nhset*2);
  elseif isequal(mfunc, 'median_dpss')
    pow  = zeros(nchan, nfreqboi*ntaper, nhset);
  else
    pow  = zeros(nchan,nfreqboi,nhset);
  end
  for ih = 1:nhset % loop resampling ratios hset
    upow = zeros(nchan,nfreqboi);
    dpow = zeros(nchan,nfreqboi);
    [n, d] = rat(hset(ih)); % n > d
    for k = 0 : nwindow-1 % loop #subset
      % pair-wised resampling
      subset_start = subset_dist*k + 1;
      subset_end = subset_start+windownsample-1;
      subdat = dat(:,subset_start:subset_end);
      udat = resample(subdat', n, d)'; % upsample
      ddat = resample(subdat', d, n)'; % downsample

      % compute auto-power of resampled data
      % upsampled
      postpad = ceil(round(pad * fsample) - size(udat,2));
      if postpad<0
        ft_error('the requested amount of zero-padding is < 0, this is not possible');
      end
      siz = [size(dat,1) endnsample ntaper];
      tmp = complex(zeros(siz), zeros(siz));
      for itap = 1:ntaper
        % fft of zero-padded tapered data segment
        tmp(:,:,itap) = fft(ft_preproc_padding(bsxfun(@times,udat,tap{ih,1}(itap,:)), padtype, 0, postpad), endnsample, 2);
      end
      if ntaper>1 && ~isequal(mfunc, 'median_dpss')
        tmp = adaptspec_dpss(tmp, lambda{ih,1}, adaptflag);
        tmp = tmp(:,freqboi) .* (2 ./ endnsample);
      else
        tmp = tmp(:,freqboi,:);
        tmp = tmp .* sqrt(2 ./ endnsample);
        tmp = abs(tmp).^2;
      end
      upow = upow + tmp;
      % downsampled
      postpad = ceil(round(pad * fsample) - size(ddat,2));
      tmp = complex(zeros(siz), zeros(siz));
      for itap = 1:ntaper
        % fft of zero-padded tapered data segment
        tmp(:,:,itap) = fft(ft_preproc_padding(bsxfun(@times,ddat,tap{ih+nhset,1}(itap,:)), padtype, 0, postpad), endnsample, 2);
      end
      if ntaper>1 && ~isequal(mfunc, 'median_dpss')
        tmp = adaptspec_dpss(tmp, lambda{ih+nhset,1}, adaptflag);
        tmp = tmp(:,freqboi) .* (2 ./ endnsample);
      else
        tmp = tmp(:,freqboi,:);
        tmp = tmp .* sqrt(2 ./ endnsample);
        tmp = abs(tmp).^2;
      end
      dpow = dpow + tmp;
    end % loop #subset

    % average across subsets for noise filtering
    upow = upow / nwindow;
    dpow = dpow / nwindow;

    if isequal(mfunc, 'none_individual') || isequal(mfunc, 'median_nomean')
      % keep the individual stretched/compressed spectra
      pow(:,:,ih)       = upow;
      pow(:,:,ih+nhset) = dpow;
    else
      % meanetric mean for relocating oscillatory component
      pow(:,:,ih) = sqrt(upow(:,:).*dpow(:,:));
    end
  end % loop resampling ratios hset

  switch mfunc
    case 'median'
      % median across resampling factors, discard the frequencies that are beyond the bandwidth after upsampling
      freqnan = freqoi(:)*hset>fsample/2;
      pow(repmat(shiftdim(freqnan, -1), [size(pow,1) 1 1])) = nan;  
      spectrum = median(pow,3,'omitnan');
    case 'median_nomean'
      % this has been suggested by Andrew Quinn in his talk at cuttingEEG,
      % with 'simple' spectra it may be effective in reducing the sidelobe
      % leakage in the residuals. Note JM: it comes at the price that the
      % cancellation of the scaling operation the stretch/compress (with
      % h and 1/h may be less effective.
      freqnan = freqoi(:)*[hset 1./hset]>fsample/2;
      pow(repmat(shiftdim(freqnan, -1), [size(pow,1) 1 1])) = nan;  
      spectrum = median(pow,3,'omitnan');
    case 'median_dpss'
      % try this: do a median for each Slepian, and then average across tapers
      freqnan = freqoi(:)*hset>fsample/2;
      freqnan = repmat(freqnan, [ntaper 1]);
      pow(repmat(shiftdim(freqnan, -1), [size(pow,1) 1 1])) = nan;  
      pow = reshape(pow, [nchan nfreqoi ntaper nhset]);
      spectrum = mean(median(pow,4,'omitnan'),3);
    case 'trimmean'
      spectrum = trimmean(pow, 0.1, 3);
    case {'none_mean' 'none_individual'}
      spectrum = permute(pow, [3 1 2]);
      hset     = [hset 1./hset];
  end


elseif strcmp(output,'original')
  %%%%%%%% FFT %%%%%%%%%%

  pow  = zeros(nchan,nfreqboi);
  for k = 0 : nwindow-1 % loop #subset
    % pair-wised resampling
    subset_start = subset_dist*k + 1;
    subset_end = subset_start+windownsample-1;
    subdat = dat(:,subset_start:subset_end);
    % compute auto-power
    postpad = ceil(round(pad * fsample) - size(subdat,2));

    siz = [size(dat,1) endnsample ntaper];
    tmp = complex(zeros(siz), zeros(siz));
    for itap = 1:ntaper
      % fft of zero-padded tapered data segment
      tmp(:,:,itap) = fft(ft_preproc_padding(bsxfun(@times,subdat,tap{1}(itap,:)), padtype, 0, postpad), endnsample, 2);
    end
    if ntaper>1
      tmp = adaptspec_dpss(tmp, lambda{ih+nhset,1}, adaptflag);
      tmp = tmp(:,freqboi) .* (2 ./ endnsample);
    else
      tmp = tmp(:,freqboi);
      tmp = tmp .* sqrt(2 ./ endnsample);
      tmp = abs(tmp).^2;
    end
    pow = pow + tmp;
  end % loop #subset

  % average across subsets for noise filtering
  spectrum = pow / nwindow;

end % if output

%spectrum = reshape(vertcat(spectrum{:}),[nchan ntaper nfreqboi]); % collecting in a cell-array and later reshaping provides significant speedups
if ndims(spectrum)==2
  spectrum = shiftdim(spectrum, -1);
  ntaper   = ones(1, numel(freqboi));
else
  ntaper   = ones(1, numel(freqboi)).*size(spectrum,1);
end

% remember the current input arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
previous_argin = current_argin;
previous_tap   = tap;
previous_lambda = lambda;
