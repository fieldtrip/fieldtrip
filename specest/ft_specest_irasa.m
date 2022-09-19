function [spectrum, ntaper, freqoi] = ft_specest_irasa(dat, time, varargin)

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
persistent previous_argin previous_tap

% get the optional input arguments
freqoi    = ft_getopt(varargin, 'freqoi', 'all');
output    = ft_getopt(varargin, 'output', 'fractal');
pad       = ft_getopt(varargin, 'pad');
padtype   = ft_getopt(varargin, 'padtype', 'zero');
polyorder = ft_getopt(varargin, 'polyorder', 1);
verbose   = ft_getopt(varargin, 'verbose', true);
fbopt     = ft_getopt(varargin, 'feedback');

% check output option
if ~strcmp(output, {'fractal','original'})
  ft_error('The current version ft_specest_irasa outputs ''fractal'' or ''original'' power only. For more information about the update, see https://www.fieldtriptoolbox.org/example/irasa/');
end

% this does not work on integer data
dat = cast(dat, 'double');
if ~isa(dat, 'double') && ~isa(dat, 'single')
  dat = cast(dat, 'double');
end

% dat param
[nchan,ndatsample] = size(dat);
fsample = 1./mean(diff(time));

% subset param
subset_nsample = 2^floor(log2(ndatsample*0.9)); % the number of sub-segments in samples is the power of 2 that does not exceed 90% of input data
subset_num     = 10;                            % the number of sub-segments
subset_dist    = floor((ndatsample  - subset_nsample)/(subset_num - 1)); % distance between sub-segements in sample, to evenly distribute the sub-subsegments within the total length of input data

% resampling ratio
hset  = 1.1:0.05:1.9;
nhset = length(hset);

% remove polynomial fit from the data -> default is demeaning
if polyorder >= 0
  dat = ft_preproc_polyremoval(dat, polyorder, 1, ndatsample);
end

% zero padding
if isempty(pad)
  % if no padding is specified this is set equal to the current data length
  pad = subset_nsample/fsample;
end
if round(pad * fsample) < subset_nsample
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

% determine whether tapers need to be recomputed
current_argin = {output, time, endnsample, freqoi}; % reasoning: if time and endnsample are equal, it's the same length trial, if the rest is equal then the requested output is equal
if isequal(current_argin, previous_argin)
  % don't recompute tapers
  tap = previous_tap;
else
  if strcmp(output, 'fractal')
    % (re)compute tapers, 1:mid are upsample tapers, mid+1:end are downsample tapers
    tap = cell(nhset*2,1);
    for ih = 1:nhset
      [n, d] = rat(hset(ih)); % n > d
      tmp = hanning(size(resample(zeros(subset_nsample,nchan), n, d),1))';
      tap{ih,1} = tmp./norm(tmp, 'fro');% for upsampled subsets
      tmp = hanning(size(resample(zeros(subset_nsample,nchan), d, n),1))';
      tap{ih+nhset,1} = tmp./norm(tmp, 'fro');% for downsampled subsets
    end
  elseif strcmp(output, 'original')
    tap = hanning(subset_nsample)';
    tap = tap./norm(tap, 'fro');
  end
end

% set ntaper
if strcmp(output, 'fractal')
  ntaper = repmat(size(tap,2),nfreqoi,1); % pretend there's only one taper
elseif strcmp(output, 'original')
  ntaper = repmat(size(tap,1),nfreqoi,1);
end

% feedback of computing progress
if isempty(fbopt)
  fbopt.i = 1;
  fbopt.n = 1;
end
str = sprintf('nfft: %d samples, datalength: %d samples, %d tapers',endnsample,ndatsample,ntaper(1));
st  = dbstack;
if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis')
  % specest_mtmfft has been called by ft_freqanalysis, meaning that ft_progress has been initialised
  ft_progress(fbopt.i./fbopt.n, ['processing trial %d/%d ',str,'\n'], fbopt.i, fbopt.n);
elseif verbose
  fprintf([str, '\n']);
end

% compute irasa or fft
spectrum = cell(ntaper(1),1);
for itap = 1:ntaper(1)
  %%%%%%%% IRASA %%%%%%%%%%
  if strcmp(output,'fractal')
    pow  = zeros(nchan,nfreqboi,nhset);
    for ih = 1:nhset % loop resampling ratios hset
      upow = zeros(nchan,nfreqboi);
      dpow = zeros(nchan,nfreqboi);
      [n, d] = rat(hset(ih)); % n > d
      for k = 0 : subset_num-1 % loop #subset
        % pair-wised resampling
        subset_start = subset_dist*k + 1;
        subset_end = subset_start+subset_nsample-1;
        subdat = dat(:,subset_start:subset_end);
        udat = resample(subdat', n, d)'; % upsample
        ddat = resample(subdat', d, n)'; % downsample
        
        % compute auto-power of resampled data
        % upsampled
        postpad = ceil(round(pad * fsample) - size(udat,2));
        tmp = fft(ft_preproc_padding(bsxfun(@times,udat,tap{ih,1}), padtype, 0, postpad), endnsample, 2);
        tmp = tmp(:,freqboi);
        tmp = tmp .* sqrt(2 ./ endnsample);
        tmp = abs(tmp).^2;
        upow = upow + tmp;
        % downsampled
        postpad = ceil(round(pad * fsample) - size(ddat,2));
        tmp = fft(ft_preproc_padding(bsxfun(@times,ddat,tap{ih+nhset,1}), padtype, 0, postpad), endnsample, 2);
        tmp = tmp(:,freqboi);
        tmp = tmp .* sqrt(2 ./ endnsample);
        tmp = abs(tmp).^2;
        dpow = dpow + tmp;
      end % loop #subset
      
      % average across subsets for noise filtering
      upow = upow / subset_num;
      dpow = dpow / subset_num;
      
      % geometric mean for relocating oscillatory component
      pow(:,:,ih) = sqrt(upow.*dpow);
    end % loop resampling ratios hset
    
    % median across resampling factors
    spectrum{itap} = median(pow,3);
    
    %%%%%%%% FFT %%%%%%%%%%
  elseif strcmp(output,'original')
    pow  = zeros(nchan,nfreqboi);
    for k = 0 : subset_num-1 % loop #subset
      % pair-wised resampling
      subset_start = subset_dist*k + 1;
      subset_end = subset_start+subset_nsample-1;
      subdat = dat(:,subset_start:subset_end);
      % compute auto-power
      postpad = ceil(round(pad * fsample) - size(subdat,2));
      tmp = fft(ft_preproc_padding(bsxfun(@times,subdat,tap(itap,:)), padtype, 0, postpad), endnsample, 2);
      tmp = tmp(:,freqboi);
      tmp = tmp .* sqrt(2 ./ endnsample);
      tmp = abs(tmp).^2;
      pow = pow + tmp;
    end % loop #subset
    
    % average across subsets for noise filtering
    spectrum{itap} = pow / subset_num;
    
  end % if output
end % for taper

spectrum = reshape(vertcat(spectrum{:}),[nchan ntaper(1) nfreqboi]); % collecting in a cell-array and later reshaping provides significant speedups
spectrum = permute(spectrum, [2 1 3]);

% remember the current input arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
previous_argin = current_argin;
previous_tap   = tap;
