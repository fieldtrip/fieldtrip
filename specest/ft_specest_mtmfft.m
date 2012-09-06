function [spectrum,ntaper,freqoi] = ft_specest_mtmfft(dat, time, varargin)

% FT_SPECEST_MTMFFT computes a fast Fourier transform using multitapering with
% the DPSS sequence or using a variety of single tapers
%
% Use as
%   [spectrum,freqoi] = specest_mtmfft(dat,time...)
% where
%   dat      = matrix of chan*sample
%   time     = vector, containing time in seconds for each sample
%   spectrum = matrix of taper*chan*freqoi of fourier coefficients
%   ntaper   = vector containing number of tapers per element of freqoi
%   freqoi   = vector of frequencies in spectrum
%
% Optional arguments should be specified in key-value pairs and can include:
%   taper      = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%   pad        = number, total length of data after zero padding (in seconds)
%   freqoi     = vector, containing frequencies of interest
%   tapsmofrq  = the amount of spectral smoothing through multi-tapering. Note: 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box
%   dimord     = 'tap_chan_freq_time' (default) or 'chan_time_freqtap' for memory efficiency (only when use variable number slepian tapers)
%   polyorder  = number, the order of the polynomial to fitted to and removed from the data
%                  prior to the fourier transform (default = 0 -> remove DC-component)
%   taperopt   = additional taper options to be used in the WINDOW function, see WINDOW
%   verbose    = output progress to console (0 or 1, default 1)
%
%
% See also FT_FREQANALYSIS, FT_SPECEST_MTMCONVOL, FT_SPECEST_TFR, FT_SPECEST_HILBERT, FT_SPECEST_WAVELET

% Copyright (C) 2010, Donders Institute for Brain, Cognition and Behaviour
%
% $Id$

% these are for speeding up computation of tapers on subsequent calls
persistent previous_argin previous_tap


% get the optional input arguments
taper     = ft_getopt(varargin, 'taper'); if isempty(taper), error('You must specify a taper'); end
pad       = ft_getopt(varargin, 'pad');
freqoi    = ft_getopt(varargin, 'freqoi', 'all');
tapsmofrq = ft_getopt(varargin, 'tapsmofrq');
dimord    = ft_getopt(varargin, 'dimord', 'tap_chan_freq_time');
fbopt     = ft_getopt(varargin, 'feedback');
verbose   = ft_getopt(varargin, 'verbose', true);
polyorder = ft_getopt(varargin, 'polyorder', 0);
tapopt    = ft_getopt(varargin, 'taperopt');


if isempty(fbopt),
  fbopt.i = 1;
  fbopt.n = 1;
end

% throw errors for required input
if isempty(tapsmofrq) && strcmp(taper, 'dpss')
  error('you need to specify tapsmofrq when using dpss tapers')
end

% Set n's
[nchan,ndatsample] = size(dat);

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
postpad    = zeros(1,ceil((pad - dattime) * fsample));
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;                   % total time in seconds of padded data



% Set freqboi and freqoi
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
if numel(tapsmofrq)~=1 && (numel(tapsmofrq)~=nfreqoi)
  error('tapsmofrq needs to contain a smoothing parameter for every frequency when requesting variable number of slepian tapers')
end



% determine whether tapers need to be recomputed
current_argin = {time, postpad, taper, tapsmofrq, freqoi, tapopt}; % reasoning: if time and postpad are equal, it's the same length trial, if the rest is equal then the requested output is equal
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
          warning_once('using only one taper for specified smoothing');
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
      tap = sine_taper(ndatsample, ndatsample*(tapsmofrq./fsample))';
      tap = tap(1:(end-1), :); % remove the last taper
      
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
if ~(strcmp(taper,'dpss') && numel(tapsmofrq)>1) % variable number of slepian tapers not requested
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
    angletransform(ifreqoi) = angle(complex(coswav,sinwav));
  end
  angletransform = repmat(angletransform,[nchan,1]);
end




% compute fft
if ~(strcmp(taper,'dpss') && numel(tapsmofrq)>1) % ariable number of slepian tapers not requested
  str = sprintf('nfft: %d samples, datalength: %d samples, %d tapers',endnsample,ndatsample,ntaper(1));
  [st, cws] = dbstack;
  if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis')
    % specest_mtmfft has been called by ft_freqanalysis, meaning that ft_progress has been initialised
    ft_progress(fbopt.i./fbopt.n, ['processing trial %d/%d ',str,'\n'], fbopt.i, fbopt.n);
  else
    fprintf([str, '\n']);
  end
  spectrum = cell(ntaper(1),1);
  for itap = 1:ntaper(1)
    dum = transpose(fft(transpose([dat .* repmat(tap(itap,:),[nchan, 1]) repmat(postpad,[nchan, 1])]))); % double explicit transpose to speedup fft
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
      spectrum = complex(zeros([max(ntaper) nchan nfreqoi]));
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
          dum = transpose(fft(transpose([dat .* repmat(tap{ifreqoi}(itap,:),[nchan, 1]) repmat(postpad,[nchan, 1])]))); % double explicit transpose to speedup fft
          dum = dum(:,freqboi(ifreqoi));
          % phase-shift according to above angles
          if timedelay ~= 0
            dum = dum .* exp(-1i*angletransform);
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
          dum = transpose(fft(transpose([dat .* repmat(tap{ifreqoi}(itap,:),[nchan, 1]) repmat(postpad,[nchan, 1])]))); % double explicit transpose to speedup fft
          dum = dum(:,freqboi(ifreqoi));
          % phase-shift according to above angles
          if timedelay ~= 0
            dum = dum .* exp(-1i*angletransform);
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
% tapers for Matlab 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin)
tap = dpss(double(a), double(b), varargin{:});
