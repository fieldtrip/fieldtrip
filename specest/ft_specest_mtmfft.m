function [spectrum,ntaper,freqoi] = ft_specest_mtmfft(dat, time, varargin)

% FT_SPECEST_MTMFFT computes a fast Fourier transform using multitapering with
% multiple tapers from the DPSS sequence or using a variety of single tapers.
%
% Use as
%   [spectrum,ntaper,freqoi] = ft_specest_mtmfft(dat,time...)
% where
%   dat        = matrix of chan*sample
%   time       = vector, containing time in seconds for each sample
%   spectrum   = matrix of taper*chan*freqoi of fourier coefficients
%   ntaper     = vector containing number of tapers per element of freqoi
%   freqoi     = vector of frequencies in spectrum
%
% Optional arguments should be specified in key-value pairs and can include
%   taper      = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%   pad        = number, total length of data after zero padding (in seconds)
%   padtype    = string, indicating type of padding to be used (see ft_preproc_padding, default: zero)
%   freqoi     = vector, containing frequencies of interest
%   tapsmofrq  = the amount of spectral smoothing through multi-tapering. Note: 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box
%   dimord     = 'tap_chan_freq' (default) or 'chan_time_freqtap' for memory efficiency (only when use variable number slepian tapers)
%   polyorder  = number, the order of the polynomial to fitted to and removed from the data prior to the fourier transform (default = 0 -> remove DC-component)
%   taperopt   = additional taper options to be used in the WINDOW function, see WINDOW
%   verbose    = output progress to console (0 or 1, default 1)
%
% See also FT_FREQANALYSIS, FT_SPECEST_MTMCONVOL, FT_SPECEST_TFR, FT_SPECEST_HILBERT, FT_SPECEST_WAVELET

% Copyright (C) 2010, Donders Institute for Brain, Cognition and Behaviour
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
taper     = ft_getopt(varargin, 'taper'); if isempty(taper), error('You must specify a taper'); end
pad       = ft_getopt(varargin, 'pad');
padtype   = ft_getopt(varargin, 'padtype', 'zero');
freqoi    = ft_getopt(varargin, 'freqoi', 'all');
tapsmofrq = ft_getopt(varargin, 'tapsmofrq');
dimord    = ft_getopt(varargin, 'dimord', 'tap_chan_freq');
fbopt     = ft_getopt(varargin, 'feedback');
verbose   = ft_getopt(varargin, 'verbose', true);
polyorder = ft_getopt(varargin, 'polyorder', 0);
tapopt    = ft_getopt(varargin, 'taperopt');

if isempty(fbopt)
  fbopt.i = 1;
  fbopt.n = 1;
end

% throw errors for required input
if isempty(tapsmofrq) && (strcmp(taper, 'dpss') || strcmp(taper, 'sine'))
  error('you need to specify tapsmofrq when using dpss or sine tapers')
end

% this does not work on integer data
dat = cast(dat, 'double');

% Set n's
[nchan,ndatsample] = size(dat);

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
  error('the padding that you specified is shorter than the data');
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
  error('tapsmofrq needs to contain a smoothing parameter for every frequency when requesting variable number of slepian tapers')
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
  % recompute tapers
  switch taper
    
    case 'dpss'
      if numel(tapsmofrq)==1
        % create a sequence of DPSS tapers, ensure that the input arguments are double precision
        tap = double_dpss(ndatsample,ndatsample*(tapsmofrq./fsample))';
        % remove the last taper because the last slepian taper is always messy
        tap = tap(1:(end-1), :);
        
        % give error/warning about number of tapers
        if isempty(tap)
          error('datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',ndatsample/fsample,tapsmofrq,fsample/ndatsample);
        elseif size(tap,1) == 1
          ft_warning('using only one taper for specified smoothing');
        end
      elseif numel(tapsmofrq)>1
        tap = cell(1,nfreqoi);
        for ifreqoi = 1:nfreqoi
          % create a sequence of DPSS tapers, ensure that the input arguments are double precision
          currtap = double_dpss(ndatsample, ndatsample .* (tapsmofrq(ifreqoi) ./ fsample))';
          % remove the last taper because the last slepian taper is always messy
          currtap = currtap(1:(end-1), :);
          
          % give error/warning about number of tapers
          if isempty(currtap)
            error('%.3f Hz: datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',freqoi(ifreqoi), ndatsample/fsample,tapsmofrq(ifreqoi),fsample/ndatsample(ifreqoi));
          elseif size(currtap,1) == 1
            disp([num2str(freqoi(ifreqoi)) ' Hz: WARNING: using only one taper for specified smoothing'])
          end
          tap{ifreqoi} = currtap;
        end
      end
      
    case 'sine'
      if numel(tapsmofrq)==1
        % create a sequence of sine tapers, 
        tap = sine_taper(ndatsample, ndatsample*(tapsmofrq./fsample))';
        % remove the last taper 
        tap = tap(1:(end-1), :);
        
        % give error/warning about number of tapers
        if isempty(tap)
          error('datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',ndatsample/fsample,tapsmofrq,fsample/ndatsample);
        elseif size(tap,1) == 1
          ft_warning('using only one taper for specified smoothing');
        end
      elseif numel(tapsmofrq)>1
        tap = cell(1,nfreqoi);
        for ifreqoi = 1:nfreqoi
          % create a sequence of sine tapers
          currtap = sine_taper(ndatsample, ndatsample .* (tapsmofrq(ifreqoi) ./ fsample))';
          % remove the last taper because the last slepian taper is always messy
          currtap = currtap(1:(end-1), :);
          
          % give error/warning about number of tapers
          if isempty(currtap)
            error('%.3f Hz: datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',freqoi(ifreqoi), ndatsample/fsample,tapsmofrq(ifreqoi),fsample/ndatsample(ifreqoi));
          elseif size(currtap,1) == 1
            disp([num2str(freqoi(ifreqoi)) ' Hz: WARNING: using only one taper for specified smoothing'])
          end
          tap{ifreqoi} = currtap;
        end
      end          
      
    case 'sine_old'
      % to provide compatibility with the tapers being scaled (which was default
      % behavior prior to 29apr2011) yet this gave different magnitude of power
      % when comparing with slepian multi tapers
      tap = sine_taper_scaled(ndatsample, ndatsample*(tapsmofrq./fsample))';
      tap = tap(1:(end-1), :); % remove the last taper
      
    case 'alpha'
      error('not yet implemented');
      
    case 'hanning'
      tap = hanning(ndatsample)';
      tap = tap./norm(tap, 'fro');
      
    otherwise
      % create the taper and ensure that it is normalized
      if isempty(tapopt) % some windowing functions don't support nargin>1, and window.m doesn't check it
        tap = window(taper, ndatsample)';
      else
        tap = window(taper, ndatsample, tapopt)';
      end
      tap = tap ./ norm(tap,'fro');
      
  end % switch taper
end % isequal currargin

% set ntaper
if ~((strcmp(taper,'dpss') || strcmp(taper,'sine')) && numel(tapsmofrq)>1) % variable number of slepian tapers not requested
  ntaper = repmat(size(tap,1),nfreqoi,1);
else % variable number of slepian tapers requested
  ntaper = cellfun(@size,tap,repmat({1},[1 nfreqoi]));
end

% determine phase-shift so that for all frequencies angle(t=0) = 0
timedelay = time(1);
if timedelay ~= 0
  angletransform = complex(zeros(1,nfreqoi));
  for ifreqoi = 1:nfreqoi
    missedsamples = round(timedelay * fsample);
    % determine angle of freqoi if oscillation started at 0
    % the angle of wavelet(cos,sin) = 0 at the first point of a cycle, with sine being in upgoing flank, which is the same convention as in mtmconvol
    anglein = (missedsamples) .* ((2.*pi./fsample) .* freqoi(ifreqoi));
    coswav  = cos(anglein);
    sinwav  = sin(anglein);
    angletransform(ifreqoi) = atan2(sinwav, coswav);
  end
  angletransform = repmat(angletransform,[nchan,1]);
end

% compute fft
if ~((strcmp(taper,'dpss') || strcmp(taper,'sine')) && numel(tapsmofrq)>1) % ariable number of slepian tapers not requested
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
    dum = fft(ft_preproc_padding(bsxfun(@times,dat,tap(itap,:)), padtype, 0, postpad),[], 2);
    dum = dum(:,freqboi);
    % phase-shift according to above angles
    if timedelay ~= 0
      dum = dum .* exp(-1i*angletransform);
    end
    dum = dum .* sqrt(2 ./ endnsample);
    spectrum{itap} = dum;
  end
  
  spectrum = reshape(vertcat(spectrum{:}),[nchan ntaper(1) nfreqboi]);% collecting in a cell-array and later reshaping provides significant speedups
  spectrum = permute(spectrum, [2 1 3]);
  
  
else % variable number of slepian tapers requested
  switch dimord
    
    case 'tap_chan_freq' % default
      % start fft'ing
      spectrum = complex(NaN([max(ntaper) nchan nfreqoi]));
      for ifreqoi = 1:nfreqoi
        str = sprintf('nfft: %d samples, datalength: %d samples, frequency %d (%.2f Hz), %d tapers',endnsample,ndatsample,ifreqoi,freqoi(ifreqoi),ntaper(ifreqoi));
        [st, cws] = dbstack;
        if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis') && verbose
          % specest_mtmconvol has been called by ft_freqanalysis, meaning that ft_progress has been initialised
          ft_progress(fbopt.i./fbopt.n, ['processing trial %d, ',str,'\n'], fbopt.i);
        elseif verbose
          fprintf([str, '\n']);
        end
        for itap = 1:ntaper(ifreqoi)
          
          dum = fft(ft_preproc_padding(bsxfun(@times,dat,tap{ifreqoi}(itap,:)), padtype, 0, postpad), [], 2);
          
          dum = dum(:,freqboi(ifreqoi));
          % phase-shift according to above angles
          if timedelay ~= 0
            dum = dum .* exp(-1i*angletransform(:,ifreqoi));
          end
          dum = dum .* sqrt(2 ./ endnsample);
          spectrum(itap,:,ifreqoi) = dum;
        end
      end % for nfreqoi
      
    case 'chan_freqtap' % memory efficient representation
      % create tapfreqind
      freqtapind = cell(1,nfreqoi);
      tempntaper = [0; cumsum(ntaper(:))];
      for ifreqoi = 1:nfreqoi
        freqtapind{ifreqoi} = tempntaper(ifreqoi)+1:tempntaper(ifreqoi+1);
      end
      
      % start fft'ing
      spectrum = complex(zeros([nchan sum(ntaper)]));
      for ifreqoi = 1:nfreqoi
        str = sprintf('nfft: %d samples, datalength: %d samples, frequency %d (%.2f Hz), %d tapers',endnsample,ndatsample,ifreqoi,freqoi(ifreqoi),ntaper(ifreqoi));
        [st, cws] = dbstack;
        if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis') && verbose
          % specest_mtmconvol has been called by ft_freqanalysis, meaning that ft_progress has been initialised
          ft_progress(fbopt.i./fbopt.n, ['processing trial %d, ',str,'\n'], fbopt.i);
        elseif verbose
          fprintf([str, '\n']);
        end
        for itap = 1:ntaper(ifreqoi)
          
          dum = fft(ft_preproc_padding(bsxfun(@times,dat,tap{ifreqoi}(itap,:)), padtype, 0, postpad), [], 2);
          
          dum = dum(:,freqboi(ifreqoi));
          % phase-shift according to above angles
          if timedelay ~= 0
            dum = dum .* exp(-1i*angletransform(:,ifreqoi));
          end
          dum = dum .* sqrt(2 ./ endnsample);
          spectrum(:,freqtapind{ifreqoi}(itap)) = dum;
        end
      end % for nfreqoi
  end % switch dimord
end

% remember the current input arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
previous_argin = current_argin;
previous_tap   = tap;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for MATLAB 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin)
tap = dpss(double(a), double(b), varargin{:});
