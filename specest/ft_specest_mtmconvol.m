function [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(dat, time, varargin)

% SPECEST_MTMCONVOL performs wavelet convolution in the time domain 
% by multiplication in the frequency domain
%
% Use as
%   [spectrum,freqoi,timeoi] = specest_mtmconvol(dat,time,...)
% where
%   dat      = matrix of chan*sample
%   time     = vector, containing time in seconds for each sample
%   spectrum = matrix of ntaper*chan*freqoi*timeoi of fourier coefficients
%   ntaper   = vector containing the number of tapers per freqoi
%   freqoi   = vector of frequencies in spectrum
%   timeoi   = vector of timebins in spectrum
%
% Optional arguments should be specified in key-value pairs and can include:
%   taper     = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%   pad       = number, indicating time-length of data to be padded out to in seconds
%   timeoi    = vector, containing time points of interest (in seconds)
%   timwin    = vector, containing length of time windows (in seconds)
%   freqoi    = vector, containing frequencies (in Hz)
%   tapsmofrq = number, the amount of spectral smoothing through multi-tapering. Note: 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box
%   dimord    = 'tap_chan_freq_time' (default) or 'chan_time_freqtap' for
%                memory efficiency
%   verbose   = output progress to console (0 or 1, default 1)
%
% See also SPECEST_MTMFFT, SPECEST_CONVOL, SPECEST_HILBERT, SPECEST_WAVELET

% Copyright (C) 2010, Donders Institute for Brain, Cognition and Behaviour
%
% $Log$

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'taper','pad','timeoi','timwin','freqoi','tapsmofrq','dimord','feedback','verbose','polyremoval'});
taper     = keyval('taper',       varargin); if isempty(taper),    taper   = 'dpss';                  end
pad       = keyval('pad',         varargin);
timeoi    = keyval('timeoi',      varargin); if isempty(timeoi),   timeoi  = 'all';                   end
timwin    = keyval('timwin',      varargin);
freqoi    = keyval('freqoi',      varargin); if isempty(freqoi),   freqoi  = 'all';                   end
tapsmofrq = keyval('tapsmofrq',   varargin);
dimord    = keyval('dimord',      varargin); if isempty(dimord),   dimord  = 'tap_chan_freq_time';    end
fbopt     = keyval('feedback',    varargin);
verbose   = keyval('verbose',     varargin); if isempty(verbose),  verbose = 1;                       end
polyorder = keyval('polyremoval', varargin); if isempty(polyorder), polyorder = 1; end

if isempty(fbopt),
  fbopt.i = 1;
  fbopt.n = 1;
end

% throw errors for required input
if isempty(tapsmofrq) && strcmp(taper, 'dpss')
  error('you need to specify tapsmofrq when using dpss tapers')
end
if isempty(timwin)
  error('you need to specify timwin')
elseif (length(timwin) ~= length(freqoi) && ~strcmp(freqoi,'all'))
  error('timwin should be of equal length as freqoi')
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

% Zero padding
if round(pad * fsample) < ndatsample
  error('the padding that you specified is shorter than the data');
end
if isempty(pad) % if no padding is specified padding is equal to current data length
  pad = dattime;
end
postpad = zeros(1,round((pad - dattime) * fsample));
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;            % total time in seconds of padded data

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

% set number of samples per time-window (timwin is in seconds)
timwinsample = round(timwin .* fsample);

% Compute tapers per frequency, multiply with wavelets and compute their fft
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
        error('%.3f Hz: datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',freqoi(ifreqoi), timwinsample(ifreqoi)/fsample,tapsmofrq(ifreqoi),fsample/timwinsample(ifreqoi));
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
      tap = window(taper, timwinsample(ifreqoi))';
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
      %       end;
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

% Switch between memory efficient representation or intuitive default representation
switch dimord
        
  case 'tap_chan_freq_time' % default
    % compute fft, major speed increases are possible here, depending on which matlab is being used whether or not it helps, which mainly focuses on orientation of the to be fft'd matrix
    datspectrum = transpose(fft(transpose([dat repmat(postpad,[nchan, 1])]))); % double explicit transpose to speedup fft
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
          dum = fftshift(transpose(ifft(transpose(datspectrum .* repmat(wltspctrm{ifreqoi}(itap,:),[nchan 1])))),2); % double explicit transpose to speedup fft
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
    freqtapind = [];
    tempntaper = [0; cumsum(ntaper(:))];
    for ifreqoi = 1:nfreqoi
      freqtapind{ifreqoi} = tempntaper(ifreqoi)+1:tempntaper(ifreqoi+1);
    end
    
    % start fft'ing
    datspectrum = transpose(fft(transpose([dat repmat(postpad,[nchan, 1])]))); % double explicit transpose to speedup fft
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
        dum = fftshift(transpose(ifft(transpose(datspectrum .* repmat(wltspctrm{ifreqoi}(itap,:),[nchan 1])))),2); % double explicit transpose to speedup fft
        tmp = complex(nan(nchan,ntimeboi),nan(nchan,ntimeboi));
        tmp(:,reqtimeboiind) = dum(:,reqtimeboi);
        tmp = tmp .* sqrt(2 ./ timwinsample(ifreqoi));
        spectrum(:,:,freqtapind{ifreqoi}(itap)) = tmp;
      end
    end % for nfreqoi
end % switch dimord

% % below code does the exact same as above, but without the trick of converting to cell-arrays for speed increases. however, when there is a huge variability in number of tapers per freqoi
% % than this approach can benefit from the fact that the array can be precreated containing nans
% % compute fft, major speed increases are possible here, depending on which matlab is being used whether or not it helps, which mainly focuses on orientation of the to be fft'd matrix
% datspectrum = transpose(fft(transpose([dat repmat(postpad,[nchan, 1])]))); % double explicit transpose to speedup fft
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
% % compute fft, major speed increases are possible here, depending on which matlab is being used whether or not it helps, which mainly focuses on orientation of the to be fft'd matrix
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
%         dum = fftshift(ifft(datspectrum(ichan,:) .* wltspctrm{ifreqoi}(itap,:),[],2)); % fftshift is necessary because of post zero-padding, not necessary when pre-padding
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
% tapers for Matlab 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin)
tap = dpss(double(a), double(b), varargin{:});
