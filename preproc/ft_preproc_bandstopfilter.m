function [filt] = ft_preproc_bandstopfilter(dat,Fs,Fbp,N,type,dir,instabilityfix,df,wintype,dev,plotfiltresp,usefftfilt)

% FT_PREPROC_BANDSTOPFILTER applies a band-stop filter to the data and thereby
% removes the spectral components in the specified frequency band
%
% Use as
%   [filt] = ft_preproc_bandstopfilter(dat, Fsample, Fbp, N, type, dir, instabilityfix)
% where
%   dat        data matrix (Nchans X Ntime)
%   Fsample    sampling frequency in Hz
%   Fbp        frequency band, specified as [Fhp Flp]
%   N          optional filter order, default is 4 (but) or dependent upon
%              frequency band and data length (fir/firls)
%   type       optional filter type, can be
%                'but' Butterworth IIR filter (default)
%                'firws' windowed sinc FIR filter
%                'fir' FIR filter using MATLAB fir1 function
%                'firls' FIR filter using MATLAB firls function (requires MATLAB Signal Processing Toolbox)
%                'brickwall' Frequency-domain filter using MATLAB FFT and iFFT function
%   dir        optional filter direction, can be
%                'onepass'         forward filter only
%                'onepass-reverse' reverse filter only, i.e. backward in time
%                'twopass'         zero-phase forward and reverse filter (default except for firws)
%                'twopass-reverse' zero-phase reverse and forward filter
%                'twopass-average' average of the twopass and the twopass-reverse
%                'onepass-zerophase' zero-phase forward filter with delay compensation (default for firws, linear-phase symmetric FIR only)
%                'onepass-minphase' minimum-phase converted forward filter (non-linear!, firws only)
%   instabilityfix optional method to deal with filter instabilities
%                'no'       only detect and give error (default)
%                'reduce'   reduce the filter order
%                'split'    split the filter in two lower-order filters, apply sequentially
%   df         optional transition width (firws)
%   wintype    optional window type (firws), can be
%                'hann'                 (max passband deviation 0.0063 [0.63%], stopband attenuation -44dB)
%                'hamming' (default)    (max passband deviation 0.0022 [0.22%], stopband attenuation -53dB)
%                'blackman'             (max passband deviation 0.0002 [0.02%], stopband attenuation -74dB)
%                'kaiser'
%   dev        optional max passband deviation/stopband attenuation (firws with kaiser window, default = 0.001 [0.1%, -60 dB])
%   plotfiltresp optional, 'yes' or 'no', plot filter responses (firws, default = 'no')
%   usefftfilt optional, 'yes' or 'no', use fftfilt instead of filter (firws, default = 'no')
%
% Note that a one- or two-pass filter has consequences for the
% strength of the filter, i.e. a two-pass filter with the same filter
% order will attenuate the signal twice as strong.
%
% Further note that the filter type 'brickwall' filters in the frequency domain,
% but may have severe issues. For instance, it has the implication that the time
% domain signal is periodic. Another issue pertains to that frequencies are
% not well defined over short time intervals; particularly for low frequencies.
%
% See also PREPROC

% Copyright (c) 2003-2014, Robert Oostenveld, Arjen Stolk, Andreas Widmann
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

% determine the size of the data
[nchans, nsamples] = size(dat);

% set the default filter order later
if nargin<4 || isempty(N)
  N = [];
end

% set the default filter type
if nargin<5 || isempty(type)
  type = 'but';
end

% set the default filter direction
if nargin<6 || isempty(dir)
  if strcmp(type, 'firws')
    dir = 'onepass-zerophase';
  else
    dir = 'twopass';
  end
end

% set the default method to deal with filter instabilities
if nargin<7|| isempty(instabilityfix)
  instabilityfix = 'no';
end

% Set default transition width later
if nargin < 8 || isempty(df)
  df = [];
end

% Set default window type
if nargin < 9 || isempty(wintype)
  wintype = 'hamming';
end

% Set default passband deviation/stopband attenuation for Kaiser window
if nargin < 10 || isempty(dev)
  if strcmp(wintype, 'kaiser')
    dev = 0.001;
  else
    dev = [];
  end
end

% Set default passband deviation/stopband attenuation for Kaiser window
if nargin < 11 || isempty(plotfiltresp)
  plotfiltresp = 'no';
end

% Set default filter function
if nargin < 12 || isempty(usefftfilt)
  usefftfilt = false;
else
  % convert to boolean value
  usefftfilt = istrue(usefftfilt);
end

% Filtering does not work on integer data
if ~isa(dat, 'double') && ~isa(dat, 'single')
  dat = cast(dat, 'double');
end

% preprocessing fails on channels that contain NaN
if any(isnan(dat(:)))
  ft_warning('FieldTrip:dataContainsNaN', 'data contains NaN values');
end

% Nyquist frequency
Fn = Fs/2;

% compute filter coefficients
switch type
  case 'but'
    if isempty(N)
      N = 4;
    end
    [B, A] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn], 'stop');

  case 'firws'

    % Input arguments
    if length(Fbp) ~= 2
        error('Two cutoff frequencies required.')
    end

    % Filter order AND transition width set?
    if ~isempty(N) && ~isempty(df)
        warning('firws:dfOverridesN', 'Filter order AND transition width set - transition width setting will override filter order.')
    elseif isempty(N) && isempty(df) % Default transition width heuristic
        df = fir_df(Fbp, Fs);
    end

    % Compute filter order from transition width
    [foo, maxDf] = fir_df(Fbp, Fs); %#ok<ASGLU>
    isOrderLow = false;
    if ~isempty(df)
      if df > maxDf
        error('Transition band too wide. Maximum transition width is %.2f Hz.', maxDf)
      end
      [N, dev] = firwsord(wintype, Fs, df, dev);
    else % Check filter order otherwise
      [df, dev] = invfirwsord(wintype, Fs, N, dev);
      if df > maxDf
        nOpt = firwsord(wintype, Fs, maxDf, dev);
        warning('firws:filterOrderLow', 'Filter order too low. For better results a minimum filter order of %d is recommended. Effective cutoff frequency might deviate from requested cutoff frequency.', nOpt)
        isOrderLow = true;
      end
    end

    % Window
    if strcmp(wintype, 'kaiser')
        beta = kaiserbeta(dev);
        win = windows('kaiser', N + 1, beta);
    else
        win = windows(wintype, N + 1);
    end

    % Impulse response
    Fbp = sort(Fbp);
    B = firws(N, Fbp / Fn, 'stop', win);
    A = 1;

    % Convert to minimum phase
    if strcmp(dir, 'onepass-minphase')
      B = minphaserceps(B);
    end

    % Twopass filtering
    if strncmp(dir, 'twopass', 7)
      pbDev = (dev + 1)^2 - 1;
      sbAtt = 20 * log10(dev^2);
      order = 2 * (length(B) - 1);
      isTwopass = true;
    else
      pbDev = dev;
      sbAtt = 20 * log10(dev);
      order = length(B) - 1;
      isTwopass = false;
    end

    % Reporting
    print_once(sprintf('Bandstop filtering data: %s, order %d, %s-windowed sinc FIR\n', dir, order, wintype));
    if ~isTwopass && ~isOrderLow % Do not report shifted cutoffs
      print_once(sprintf('  cutoff (-6 dB) %g Hz and %g Hz\n', Fbp(1), Fbp(2)));
      tb = [max([Fbp(1) - df / 2 0]) Fbp(1) + df / 2 Fbp(2) - df / 2 min([Fbp(2) + df / 2 Fn])]; % Transition band edges
      print_once(sprintf('  transition width %.1f Hz, passband 0-%.1f Hz, stopband %.1f-%.1f Hz, passband %.1f-%.0f Hz\n', df, tb, Fn));
    end
    if ~isOrderLow
      print_once(sprintf('  max. passband deviation %.4f (%.2f%%), stopband attenuation %.0f dB\n', pbDev, pbDev * 100, sbAtt));
    end

    % Plot filter responses
    if strcmp(plotfiltresp, 'yes')
      plotfresp(B, [], [], Fs, dir)
    end

  case 'fir'
    if isempty(N)
      N = 3*fix(Fs / Fbp(1));
      if rem(N,2)==1,   N=N+1;    end
    end
    if N > floor( (size(dat,2) - 1) / 3)
      N=floor(size(dat,2)/3) - 2;
      if rem(N,2)==1,   N=N+1;    end
    end
    B = fir1(N, [min(Fbp)/Fn max(Fbp)/Fn], 'stop');
    A = 1;
  case 'firls' % from NUTMEG's implementation
    % Deprecated: see bug 2453
    warning('The filter type you requested is not recommended for neural signals, only proceed if you know what you are doing.')
    if isempty(N)
      N = 3*fix(Fs / Fbp(1));
    end
    if N > floor( (size(dat,2) - 1) / 3)
      N=floor(size(dat,2)/3) - 1;
    end
    f = 0:0.001:1;
    if rem(length(f),2)~=0
      f(end)=[];
    end
    z = ones(1,length(f));
    if(isfinite(min(Fbp)))
      [val,pos1] = min(abs(Fs*f/2 - min(Fbp)));
    else
      [val,pos2] = min(abs(Fs*f/2 - max(Fbp)));
      pos1=pos2;
    end
    if(isfinite(max(Fbp)))
      [val,pos2] = min(abs(Fs*f/2 - max(Fbp)));
    else
      pos2 = length(f);
    end
    z(pos1:pos2) = 0;
    A = 1;
    B = firls(N,f,z); % requires MATLAB signal processing toolbox
  case 'brickwall'
    ax = linspace(0, Fs, size(dat,2));  % frequency coefficients
    fl = nearest(ax, min(Fbp))-1;       % low cut-off frequency
    fh = nearest(ax, max(Fbp))+1;       % high cut-off frequency
    a  = 0; % suppresion rate of frequencies-not-of-interest
    f           = fft(dat,[],2);        % FFT
    f(:,fl:fh)  = a.*f(:,fl:fh);        % perform band cut-off
    filt        = 2*real(ifft(f,[],2)); % iFFT
    return
  otherwise
    error('unsupported filter type "%s"', type);
end

% demean the data before filtering
meandat = mean(dat,2);
dat = bsxfun(@minus, dat, meandat);

try
  filt = filter_with_correction(B,A,dat,dir,usefftfilt);
catch
  switch instabilityfix
    case 'no'
      rethrow(lasterror);
    case 'reduce'
      warning('backtrace', 'off')
      warning('instability detected - reducing the %dth order filter to an %dth order filter', N, N-1);
      warning('backtrace', 'on')
      filt = ft_preproc_bandstopfilter(dat,Fs,Fbp,N-1,type,dir,instabilityfix);
    case 'split'
      N1 = ceil(N/2);
      N2 = floor(N/2);
      warning('backtrace', 'off')
      warning('instability detected - splitting the %dth order filter in a sequential %dth and a %dth order filter', N, N1, N2);
      warning('backtrace', 'on')
      filt1 = ft_preproc_bandstopfilter(dat  ,Fs,Fbp,N1,type,dir,instabilityfix);
      filt  = ft_preproc_bandstopfilter(filt1,Fs,Fbp,N2,type,dir,instabilityfix);
    otherwise
      error('incorrect specification of instabilityfix');
  end % switch
end

% add the mean back to the filtered data
filt = bsxfun(@plus, filt, meandat);

