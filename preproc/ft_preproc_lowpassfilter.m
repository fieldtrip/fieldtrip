function [filt, B, A] = ft_preproc_lowpassfilter(dat, Fs, Flp, N, type, dir, instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt)

% FT_PREPROC_LOWPASSFILTER applies a low-pass filter to the data and thereby removes
% all high frequency components in the data
%
% Use as
%   [filt] = ft_preproc_lowpassfilter(dat, Fs, Flp, N, type, dir, instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt)
% where
%   dat             data matrix (Nchans X Ntime)
%   Fs              sampling frequency in Hz
%   Flp             filter frequency in Hz
%   order           optional filter order, default is 6 (but) or dependent on frequency band and data length (fir/firls)
%   type            optional filter type, can be
%                     'but'       Butterworth IIR filter (default)
%                     'firws'     FIR filter with windowed sinc
%                     'fir'       FIR filter using MATLAB fir1 function
%                     'firls'     FIR filter using MATLAB firls function (requires MATLAB Signal Processing Toolbox)
%                     'brickwall' frequency-domain filter using forward and inverse FFT
%   dir             optional filter direction, can be
%                     'onepass'                   forward filter only
%                     'onepass-reverse'           reverse filter only, i.e. backward in time
%                     'onepass-zerophase'         zero-phase forward filter with delay compensation (default for firws, linear-phase symmetric FIR only)
%                     'onepass-reverse-zerophase' zero-phase reverse filter with delay compensation
%                     'onepass-minphase'          minimum-phase converted forward filter (non-linear, only for firws)
%                     'twopass'                   zero-phase forward and reverse filter (default, except for firws)
%                     'twopass-reverse'           zero-phase reverse and forward filter
%                     'twopass-average'           average of the twopass and the twopass-reverse
%   instabilityfix  optional method to deal with filter instabilities
%                     'no'       only detect and give error (default)
%                     'reduce'   reduce the filter order
%                     'split'    split the filter in two lower-order filters, apply sequentially
%   df              optional transition width (firws)
%   wintype         optional window type (firws), can be
%                     'hamming' (default)    maximum passband deviation 0.0022 [0.22%], stopband attenuation -53dB
%                     'hann'                 maximum passband deviation 0.0063 [0.63%], stopband attenuation -44dB
%                     'blackman'             maximum passband deviation 0.0002 [0.02%], stopband attenuation -74dB
%                     'kaiser'
%   dev             optional max passband deviation/stopband attenuation (only for firws with kaiser window, default = 0.001 [0.1%, -60 dB])
%   plotfiltresp    optional, 'yes' or 'no', plot filter responses (only for firws, default = 'no')
%   usefftfilt      optional, 'yes' or 'no', use fftfilt instead of filter (only for firws, default = 'no')
%
% Note that a one- or two-pass filter has consequences for the strength of the
% filter, i.e. a two-pass filter with the same filter order will attenuate the signal
% twice as strong.
%
% Further note that the filter type 'brickwall' filters in the frequency domain,
% but may have severe issues. For instance, it has the implication that the time
% domain signal is periodic. Another issue pertains to that frequencies are
% not well defined over short time intervals; particularly for low frequencies.
%
% If the data contains NaNs, these will affect the output. With an IIR
% filter, and/or with FFT-filtering, local NaNs will spread to the whole
% time series. With a FIR filter, local NaNs will spread locally, depending
% on the filter order.
%
% See also PREPROC

% Copyright (c) 2003-2022, Robert Oostenveld, Arjen Stolk, Andreas Widmann,
% Jan-Mathijs Schoffelen
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

% demean the data before filtering
meandat = nanmean(dat,2);
dat = bsxfun(@minus, dat, meandat);

% compute filter coefficients
switch type
  case 'but'
    if isempty(N)
      N = 6;
    end
    [B, A] = butter(N, max(Flp)/Fn);
    
  case 'firws'
    % Input arguments
    if length(Flp) ~= 1
      ft_error('One cutoff frequency required.')
    end
    
    % Filter order AND transition width set?
    if ~isempty(N) && ~isempty(df)
      ft_warning('firws:dfOverridesN', 'Filter order AND transition width set - transition width setting will override filter order.')
    elseif isempty(N) && isempty(df) % Default transition width heuristic
      df = fir_df(Flp, Fs);
    end
    
    % Compute filter order from transition width
    [foo, maxDf] = fir_df(Flp, Fs); %#ok<ASGLU>
    isOrderLow = false;
    if ~isempty(df)
      if df > maxDf
        ft_error('Transition band too wide. Maximum transition width is %.2f Hz.', maxDf)
      end
      [N, dev] = firwsord(wintype, Fs, df, dev);
    else % Check filter order otherwise
      [df, dev] = invfirwsord(wintype, Fs, N, dev);
      if df > maxDf
        nOpt = firwsord(wintype, Fs, maxDf, dev);
        ft_warning('firws:filterOrderLow', 'Filter order too low. For better results a minimum filter order of %d is recommended. Effective cutoff frequency might deviate from requested cutoff frequency.', nOpt)
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
    B = firws(N, Flp / Fn, 'low', win);
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
    ft_info once
    ft_info('Lowpass filtering data: %s, order %d, %s-windowed sinc FIR\n', dir, order, wintype);
    if ~isTwopass && ~isOrderLow % Do not report shifted cutoffs
      ft_info('  cutoff (-6 dB) %g Hz\n', Flp);
      tb = [max([Flp - df / 2 0]), min([Flp + df / 2 Fn])]; % Transition band edges
      ft_info('  transition width %.1f Hz, passband 0-%.1f Hz, stopband %.1f-%.0f Hz\n', df, tb, Fn);
    end
    if ~isOrderLow
      ft_info('  maximum passband deviation %.4f (%.2f%%), stopband attenuation %.0f dB\n', pbDev, pbDev * 100, sbAtt);
    end
    
  case 'fir'
    if isempty(N)
      N = 3*fix(Fs / Flp);
    end
    if N > floor( (size(dat,2) - 1) / 3)
      N=floor(size(dat,2)/3) - 1;
    end
    B = fir1(N, max(Flp)/Fn);
    A = 1;

  case 'firls' % from NUTMEG's implementation
    % Deprecated: see bug 2453
    ft_warning('The filter type you requested is not recommended for neural signals, only proceed if you know what you are doing.')
    if isempty(N)
      N = 3*fix(Fs / Flp);
    end
    if N > floor( (size(dat,2) - 1) / 3)
      N=floor(size(dat,2)/3) - 1;
    end
    f = 0:0.001:1;
    if rem(length(f),2)~=0
      f(end)=[];
    end
    z = zeros(1,length(f));
    [val,pos1] = min(abs(Fs*f/2 - 0));
    if(isfinite(Flp))
      [val,pos2] = min(abs(Fs*f/2 - Flp));
    else
      pos2 = length(f);
    end
    z(pos1:pos2) = 1;
    A = 1;
    B = firls(N,f,z); % requires MATLAB signal processing toolbox
  case 'brickwall'
    ax = (0:(size(dat,2)-1))./(Fs/size(dat,2)); % frequency axis
    
    a     = ones(1, size(dat,2));
    fbin1 = nearest(ax, Flp);
    fbin2 = nearest(ax, Fs-Flp);
    
    a((fbin1+1):(fbin2-1)) = 0;
    
    f    = fft(dat,[],2);             % FFT
    f    = f.*a(ones(size(dat,1),1),:); % brickwall
    filt = real(ifft(f,[],2));        % iFFT
 
  otherwise
    ft_error('unsupported filter type "%s"', type);
end

% Plot filter responses
if strcmp(plotfiltresp, 'yes')
  plotfresp(B, A, [], Fs, dir)
end

if ~isequal(type, 'brickwall')
  try
    filt = filter_with_correction(B,A,dat,dir,usefftfilt);
  catch
    switch instabilityfix
      case 'no'
        rethrow(lasterror);
      case 'reduce'
        ft_warning('off','backtrace');
        ft_warning('instability detected - reducing the %dth order filter to an %dth order filter', N, N-1);
        ft_warning('on','backtrace');
        filt = ft_preproc_lowpassfilter(dat,Fs,Flp,N-1,type,dir,instabilityfix);
      case 'split'
        N1 = ceil(N/2);
        N2 = floor(N/2);
        ft_warning('off','backtrace');
        ft_warning('instability detected - splitting the %dth order filter in a sequential %dth and a %dth order filter', N, N1, N2);
        ft_warning('on','backtrace');
        filt = ft_preproc_lowpassfilter(dat ,Fs,Flp,N1,type,dir,instabilityfix);
        filt = ft_preproc_lowpassfilter(filt,Fs,Flp,N2,type,dir,instabilityfix);
      otherwise
        ft_error('incorrect specification of instabilityfix');
    end % switch
  end
end

% add the mean back to the filtered data
filt = bsxfun(@plus, filt, meandat);
