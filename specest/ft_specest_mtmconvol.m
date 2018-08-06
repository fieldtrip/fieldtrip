function [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(dat, time, varargin)

% FT_SPECEST_MTMCONVOL performs wavelet convolution in the time domain by
% multiplication in the frequency domain.
%
% Use as
%   [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(dat,time,...)
% where input
%   dat       = matrix of chan*sample
%   time      = vector, containing time in seconds for each sample
% and output
%   spectrum  = matrix of ntaper*chan*freqoi*timeoi of fourier coefficients
%   ntaper    = vector containing the number of tapers per freqoi
%   freqoi    = vector of frequencies in spectrum
%   timeoi    = vector of timebins in spectrum
%
% Optional arguments should be specified in key-value pairs and can include
%   taper     = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%   pad       = number, indicating time-length of data to be padded out to in seconds
%   padtype   = string, indicating type of padding to be used (see ft_preproc_padding, default: zero)
%   timeoi    = vector, containing time points of interest (in seconds)
%   timwin    = vector, containing length of time windows (in seconds)
%   freqoi    = vector, containing frequencies (in Hz)
%   tapsmofrq = number, the amount of spectral smoothing through multi-tapering. Note: 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box
%   dimord    = 'tap_chan_freq_time' (default) or 'chan_time_freqtap' for memory efficiency
%   verbose   = output progress to console (0 or 1, default 1)
%   taperopt  = additional taper options to be used in the WINDOW function, see WINDOW
%   polyorder = number, the order of the polynomial to fitted to and removed from the data prior to the fourier transform (default = 0 -> remove DC-component)
%
% See also FT_FREQANALYSIS, FT_SPECEST_MTMFFT, FT_SPECEST_TFR, FT_SPECEST_HILBERT, FT_SPECEST_WAVELET

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
persistent previous_argin previous_wltspctrm

% get the optional input arguments
taper     = ft_getopt(varargin, 'taper', 'dpss');
pad       = ft_getopt(varargin, 'pad');
padtype   = ft_getopt(varargin, 'padtype', 'zero');
timeoi    = ft_getopt(varargin, 'timeoi', 'all');
timwin    = ft_getopt(varargin, 'timwin');
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
  ft_error('you need to specify tapsmofrq when using dpss tapers')
end
if isempty(timwin)
  ft_error('you need to specify timwin')
elseif (length(timwin) ~= length(freqoi) && ~strcmp(freqoi,'all'))
  ft_error('timwin should be of equal length as freqoi')
end

% Set n's
[nchan, ndatsample] = size(dat);

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
postpad    = round((pad - dattime) * fsample);
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;            % total time in seconds of padded data

% Set freqboi and freqoi
freqoiinput = freqoi;
if isnumeric(freqoi) % if input is a vector
  freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1; % is equivalent to: round(freqoi .* endtime) + 1;
  freqboi   = unique(freqboi);
  freqoi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output
elseif strcmp(freqoi,'all')
  freqboilim = round([0 fsample/2] ./ (fsample ./ endnsample)) + 1;
  freqboi    = freqboilim(1):1:freqboilim(2);
  freqoi     = (freqboi-1) ./ endtime;
end
% check for freqoi = 0 and remove it, there is no wavelet for freqoi = 0
if freqoi(1)==0
  freqoi(1)  = [];
  freqboi(1) = [];
  if length(timwin) == (length(freqoi) + 1)
    timwin(1) = [];
  end
end
nfreqboi = length(freqboi);
nfreqoi  = length(freqoi);

% throw a warning if input freqoi is different from output freqoi
if isnumeric(freqoiinput)
  % check whether padding is appropriate for the requested frequency resolution
  rayl = 1/endtime;
  if any(rem(freqoiinput,rayl))
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

% Set timeboi and timeoi
timeoiinput = timeoi;
offset = round(time(1)*fsample);
if isnumeric(timeoi) % if input is a vector
  timeoi   = unique(round(timeoi .* fsample) ./ fsample);
  timeboi  = round(timeoi .* fsample - offset) + 1;
  ntimeboi = length(timeboi);
elseif strcmp(timeoi,'all') % if input was 'all'
  timeboi  = 1:length(time);
  ntimeboi = length(timeboi);
  timeoi   = time;
end

% throw a warning if input timeoi is different from output timeoi
if isnumeric(timeoiinput)
  if numel(timeoiinput) ~= numel(timeoi) % timeoi will not contain double time-bins when requested
    ft_warning('output time-bins are different from input time-bins, multiples of the same bin were requested but not given');
  else
    if any(abs(timeoiinput-timeoi) >= eps*1e6) 
      ft_warning('output time-bins are different from input time-bins');
    end
  end
end

% set number of samples per time-window (timwin is in seconds)
if numel(timwin)==1 && nfreqoi~=1
  timwin = repmat(timwin,[1 nfreqoi]);
end
timwinsample = round(timwin .* fsample);


% determine whether tapers need to be recomputed
current_argin = {time, postpad, taper, timwinsample, tapsmofrq, freqoi, timeoi, tapopt}; % reasoning: if time and postpad are equal, it's the same length trial, if the rest is equal then the requested output is equal
if isequal(current_argin, previous_argin)
  % don't recompute tapers
  wltspctrm = previous_wltspctrm;
  ntaper    = cellfun(@size,wltspctrm,repmat({1},[1 nfreqoi])');
else
  % recompute tapers per frequency, multiply with wavelets and compute their fft
  wltspctrm = cell(nfreqoi,1);
  ntaper    = zeros(nfreqoi,1);
  for ifreqoi = 1:nfreqoi
    switch taper
      case 'dpss'
        % create a sequence of DPSS tapers, ensure that the input arguments are double precision
        tap = double_dpss(timwinsample(ifreqoi), timwinsample(ifreqoi) .* (tapsmofrq(ifreqoi) ./ fsample))';
        % remove the last taper because the last slepian taper is always messy
        tap = tap(1:(end-1), :);
        
        % give error/warning about number of tapers
        if isempty(tap)
          ft_error('%.3f Hz: datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',freqoi(ifreqoi), timwinsample(ifreqoi)/fsample,tapsmofrq(ifreqoi),fsample/timwinsample(ifreqoi));
        elseif size(tap,1) == 1
          disp([num2str(freqoi(ifreqoi)) ' Hz: WARNING: using only one taper for specified smoothing'])
        end
        
      case 'sine'
        % create and remove the last taper
        tap = sine_taper(timwinsample(ifreqoi), timwinsample(ifreqoi) .* (tapsmofrq(ifreqoi) ./ fsample))';
        tap = tap(1:(end-1), :);
        
      case 'sine_old'
        % to provide compatibility with the tapers being scaled (which was default
        % behavior prior to 29apr2011) yet this gave different magnitude of power
        % when comparing with slepian multi tapers
        tap = sine_taper_scaled(timwinsample(ifreqoi), timwinsample(ifreqoi) .* (tapsmofrq(ifreqoi) ./ fsample))';
        tap = tap(1:(end-1), :); % remove the last taper
        
      case 'alpha'
        tap = alpha_taper(timwinsample(ifreqoi), freqoi(ifreqoi)./ fsample)';
        tap = tap./norm(tap)';
        
      case 'hanning'
        tap = hanning(timwinsample(ifreqoi))';
        tap = tap./norm(tap, 'fro');
        
      otherwise
        % create a single taper according to the window specification as a replacement for the DPSS (Slepian) sequence
        if isempty(tapopt) % some windowing functions don't support nargin>1, and window.m doesn't check it
          tap = window(taper, timwinsample(ifreqoi))';
        else
          tap = window(taper, timwinsample(ifreqoi),tapopt)';
        end
        tap = tap ./ norm(tap,'fro'); % make it explicit that the frobenius norm is being used
    end
    
    % set number of tapers
    ntaper(ifreqoi) = size(tap,1);
    
    % Wavelet construction
    tappad   = ceil(endnsample ./ 2) - floor(timwinsample(ifreqoi) ./ 2);
    prezero  = zeros(1,tappad);
    postzero = zeros(1,round(endnsample) - ((tappad-1) + timwinsample(ifreqoi))-1);
    
    % phase consistency: cos must always be 1  and sin must always be centered in upgoing flank, so the centre of the wavelet (untapered) has angle = 0
    anglein  = (-(timwinsample(ifreqoi)-1)/2 : (timwinsample(ifreqoi)-1)/2)'   .*  ((2.*pi./fsample) .* freqoi(ifreqoi));
    wltspctrm{ifreqoi} = complex(zeros(size(tap,1),round(endnsample)));
    
    for itap = 1:ntaper(ifreqoi)
      try % this try loop tries to fit the wavelet into wltspctrm, when its length is smaller than ndatsample, the rest is 'filled' with zeros because of above code
        % if a wavelet is longer than ndatsample, it doesn't fit and it is kept at zeros, which is translated to NaN's in the output
        % construct the complex wavelet
        coswav  = horzcat(prezero, tap(itap,:) .* cos(anglein)', postzero);
        sinwav  = horzcat(prezero, tap(itap,:) .* sin(anglein)', postzero);
        wavelet = complex(coswav, sinwav);
        % store the fft of the complex wavelet
        wltspctrm{ifreqoi}(itap,:) = fft(wavelet,[],2);
        
        %       % debug plotting
        %       figure('name',['taper #' num2str(itap) ' @ ' num2str(freqoi(ifreqoi)) 'Hz' ],'NumberTitle','off');
        %       subplot(2,1,1);
        %       hold on;
        %       plot(real(wavelet));
        %       plot(imag(wavelet),'color','r');
        %       legend('real','imag');
        %       tline = length(wavelet)/2;
        %       if mod(tline,2)==0
        %         line([tline tline],[-max(abs(wavelet)) max(abs(wavelet))],'color','g','linestyle','--')
        %       else
        %         line([ceil(tline) ceil(tline)],[-max(abs(wavelet)) max(abs(wavelet))],'color','g','linestyle','--');
        %         line([floor(tline) floor(tline)],[-max(abs(wavelet)) max(abs(wavelet))],'color','g','linestyle','--');
        %       end
        %       subplot(2,1,2);
        %       plot(angle(wavelet),'color','g');
        %       if mod(tline,2)==0,
        %         line([tline tline],[-pi pi],'color','r','linestyle','--')
        %       else
        %         line([ceil(tline) ceil(tline)],[-pi pi],'color','r','linestyle','--')
        %         line([floor(tline) floor(tline)],[-pi pi],'color','r','linestyle','--')
        %       end
      end
    end
  end
end

% Switch between memory efficient representation or intuitive default representation
switch dimord
        
  case 'tap_chan_freq_time' % default
    % compute fft, major speed increases are possible here, depending on which MATLAB is being used whether or not it helps, which mainly focuses on orientation of the to be fft'd matrix
    datspectrum = fft(ft_preproc_padding(dat, padtype, 0, postpad), [], 2);
    spectrum = cell(max(ntaper), nfreqoi);
    for ifreqoi = 1:nfreqoi
      str = sprintf('frequency %d (%.2f Hz), %d tapers', ifreqoi,freqoi(ifreqoi),ntaper(ifreqoi));
      [st, cws] = dbstack;
      if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis') && verbose
        % specest_mtmconvol has been called by ft_freqanalysis, meaning that ft_progress has been initialised
        ft_progress(fbopt.i./fbopt.n, ['trial %d, ',str,'\n'], fbopt.i);
      elseif verbose
        fprintf([str, '\n']);
      end

      for itap = 1:max(ntaper)
        % compute indices that will be used to extracted the requested fft output
        nsamplefreqoi    = timwin(ifreqoi) .* fsample;
        reqtimeboiind    = find((timeboi >=  (nsamplefreqoi ./ 2)) & (timeboi <    ndatsample - (nsamplefreqoi ./2)));
        reqtimeboi       = timeboi(reqtimeboiind);
        
        % compute datspectrum*wavelet, if there are reqtimeboi's that have data
        % create a matrix of NaNs if there is no taper for this current frequency-taper-number
        if itap > ntaper(ifreqoi)
          spectrum{itap,ifreqoi} = complex(nan(nchan,ntimeboi));
        else                                
          dum = fftshift(ifft(datspectrum .* repmat(wltspctrm{ifreqoi}(itap,:),[nchan 1]), [], 2),2); % fftshift is necessary to implement zero-phase/acyclic/acausal convolution (either here, or the wavelet should be wrapped around sample=0)
          tmp = complex(nan(nchan,ntimeboi));
          tmp(:,reqtimeboiind) = dum(:,reqtimeboi);
          tmp = tmp .* sqrt(2 ./ timwinsample(ifreqoi));
          spectrum{itap,ifreqoi} = tmp;
        end
      end
    end
    spectrum = reshape(vertcat(spectrum{:}),[nchan max(ntaper) nfreqoi ntimeboi]); % collecting in a cell-array and later reshaping provides significant speedups
    spectrum = permute(spectrum, [2 1 3 4]);
    
  case 'chan_time_freqtap' % memory efficient representation
    % create tapfreqind
    freqtapind = cell(1,nfreqoi);
    tempntaper = [0; cumsum(ntaper(:))];
    for ifreqoi = 1:nfreqoi
      freqtapind{ifreqoi} = tempntaper(ifreqoi)+1:tempntaper(ifreqoi+1);
    end
    
    % start fft'ing
    datspectrum = fft(ft_preproc_padding(dat, padtype, 0, postpad), [], 2);
    spectrum = complex(zeros([nchan ntimeboi sum(ntaper)]));
    for ifreqoi = 1:nfreqoi
      str = sprintf('frequency %d (%.2f Hz), %d tapers', ifreqoi,freqoi(ifreqoi),ntaper(ifreqoi));
      [st, cws] = dbstack;
      if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis') && verbose
        % specest_mtmconvol has been called by ft_freqanalysis, meaning that ft_progress has been initialised
        ft_progress(fbopt.i./fbopt.n, ['trial %d, ',str,'\n'], fbopt.i);
      elseif verbose
        fprintf([str, '\n']);
      end
      for itap = 1:ntaper(ifreqoi)
        % compute indices that will be used to extracted the requested fft output
        nsamplefreqoi    = timwin(ifreqoi) .* fsample;
        reqtimeboiind    = find((timeboi >=  (nsamplefreqoi ./ 2)) & (timeboi <    (ndatsample - (nsamplefreqoi ./2))));
        reqtimeboi       = timeboi(reqtimeboiind);
        
        % compute datspectrum*wavelet, if there are reqtimeboi's that have data
        dum = fftshift(ifft(datspectrum .* repmat(wltspctrm{ifreqoi}(itap,:),[nchan 1]), [], 2),2); % fftshift is necessary to implement zero-phase/acyclic/acausal convolution (either here, or the wavelet should be wrapped around sample=0)
        tmp = complex(nan(nchan,ntimeboi),nan(nchan,ntimeboi));
        tmp(:,reqtimeboiind) = dum(:,reqtimeboi);
        tmp = tmp .* sqrt(2 ./ timwinsample(ifreqoi));
        spectrum(:,:,freqtapind{ifreqoi}(itap)) = tmp;
      end
    end % for nfreqoi
end % switch dimord


% remember the current input arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
previous_argin     = current_argin;
previous_wltspctrm = wltspctrm;


% % below code does the exact same as above, but without the trick of converting to cell-arrays for speed increases. however, when there is a huge variability in number of tapers per freqoi
% % than this approach can benefit from the fact that the array can be precreated containing nans
% % compute fft, major speed increases are possible here, depending on which MATLAB is being used whether or not it helps, which mainly focuses on orientation of the to be fft'd matrix
% datspectrum = transpose(fft(transpose([dat repmat(postpad,[nchan, 1])]))); % double explicit transpose to speedup fft
% % NOTE: double explicit transpose around fft is not faster than using fft
% with dim argument (in fact, it is slower)
% spectrum = complex(nan([max(ntaper) nchan nfreqoi ntimeboi]),nan([max(ntaper) nchan nfreqoi ntimeboi])); % assumes fixed number of tapers
% for ifreqoi = 1:nfreqoi
%   fprintf('processing frequency %d (%.2f Hz), %d tapers\n', ifreqoi,freqoi(ifreqoi),ntaper(ifreqoi));
%   for itap = 1:max(ntaper)
%     % compute indices that will be used to extracted the requested fft output
%     nsamplefreqoi    = timwin(ifreqoi) .* fsample;
%     reqtimeboiind    = find((timeboi >=  (nsamplefreqoi ./ 2)) & (timeboi <    ndatsample - (nsamplefreqoi ./2)));
%     reqtimeboi       = timeboi(reqtimeboiind);
%
%     % compute datspectrum*wavelet, if there are reqtimeboi's that have data
%     % create a matrix of NaNs if there is no taper for this current frequency-taper-number
%     if itap <= ntaper(ifreqoi)
%       dum = fftshift(transpose(ifft(transpose(datspectrum .* repmat(wltspctrm{ifreqoi}(itap,:),[nchan 1])))),2); % double explicit transpose to speedup fft
%       tmp = complex(nan(nchan,ntimeboi));
%       tmp(:,reqtimeboiind) = dum(:,reqtimeboi);
%       tmp = tmp .* sqrt(2 ./ timwinsample(ifreqoi));
%       spectrum(itap,:,ifreqoi,:) = tmp;
%     else
%       break
%     end
%   end
% end




% Below the code used to implement variable amount of tapers in a different way, kept here for testing, please do not remove
% % build tapfreq vector
% tapfreq = [];
% for ifreqoi = 1:nfreqoi
%   tapfreq = [tapfreq ones(1,ntaper(ifreqoi)) * ifreqoi];
% end
% tapfreq = tapfreq(:);
%
% % compute fft, major speed increases are possible here, depending on which MATLAB is being used whether or not it helps, which mainly focuses on orientation of the to be fft'd matrix
% %spectrum = complex(nan([numel(tapfreq),nchan,ntimeboi]));
% datspectrum = fft([dat repmat(postpad,[nchan, 1])],[],2);
% spectrum = cell(numel(tapfreq), nchan, ntimeboi);
% for ifreqoi = 1:nfreqoi
%   fprintf('processing frequency %d (%.2f Hz), %d tapers\n', ifreqoi,freqoi(ifreqoi),ntaper(ifreqoi));
%   for itap = 1:ntaper(ifreqoi)
%     tapfreqind = sum(ntaper(1:ifreqoi-1)) + itap;
%     for ichan = 1:nchan
%       % compute indices that will be used to extracted the requested fft output
%       nsamplefreqoi    = timwin(ifreqoi) .* fsample;
%       reqtimeboiind    = find((timeboi >=  (nsamplefreqoi ./ 2)) & (timeboi <    ndatsample - (nsamplefreqoi ./2)));
%       reqtimeboi       = timeboi(reqtimeboiind);
%
%       % compute datspectrum*wavelet, if there are reqtimeboi's that have data
%       if ~isempty(reqtimeboi)
%         dum = fftshift(ifft(datspectrum(ichan,:) .* wltspctrm{ifreqoi}(itap,:),[],2)); 
%         %spectrum(tapfreqind,ichan,reqtimeboiind) = dum(reqtimeboi);
%         tmp = complex(nan(1,ntimeboi));
%         tmp(reqtimeboiind) = dum(reqtimeboi);
%         spectrum{tapfreqind,ichan} = tmp;
%       end
%     end
%   end
% end
% spectrum = reshape(vertcat(spectrum{:}),[numel(tapfreq) nchan ntimeboi]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for MATLAB 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin)
tap = dpss(double(a), double(b), varargin{:});
