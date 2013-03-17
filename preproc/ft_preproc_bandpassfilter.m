function [filt] = ft_preproc_bandpassfilter(dat,Fs,Fbp,N,type,dir,instabilityfix)

% FT_PREPROC_BANDPASSFILTER applies a band-pass filter to the data and thereby
% removes the spectral components in the data except for the ones in the
% specified frequency band
%
% Use as
%   [filt] = ft_preproc_bandpassfilter(dat, Fsample, Fbp, N, type, dir, instabilityfix)
% where
%   dat        data matrix (Nchans X Ntime)
%   Fsample    sampling frequency in Hz
%   Fbp        frequency band, specified as [Fhp Flp]
%   N          optional filter order, default is 4 (but) or dependent upon
%              frequency band and data length (fir/firls)
%   type       optional filter type, can be
%                'but' Butterworth IIR filter (default)
%                'fir' FIR filter using Matlab fir1 function
%                'firls' FIR filter using Matlab firls function (requires Matlab Signal Processing Toolbox)
%   dir        optional filter direction, can be
%                'onepass'         forward filter only
%                'onepass-reverse' reverse filter only, i.e. backward in time
%                'twopass'         zero-phase forward and reverse filter (default)
%                'twopass-reverse' zero-phase reverse and forward filter
%                'twopass-average' average of the twopass and the twopass-reverse
%   instabilityfix optional method to deal with filter instabilities
%                'no'       only detect and give error (default)
%                'reduce'   reduce the filter order
%                'split'    split the filter in two lower-order filters, apply sequentially
%
% Note that a one- or two-pass filter has consequences for the
% strength of the filter, i.e. a two-pass filter with the same filter
% order will attenuate the signal twice as strong.
%
% See also PREPROC

% Copyright (c) 2003-2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
if nargin<6|| isempty(dir)
  dir = 'twopass';
end

% set the default method to deal with filter instabilities
if nargin<7|| isempty(instabilityfix)
  instabilityfix = 'none';
end

% Nyquist frequency
Fn = Fs/2;

% compute filter coefficients
switch type
  case 'but'
    if isempty(N)
      N = 4;
    end
    [B, A] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn]);
  case 'fir'
    if isempty(N)
      N = 3*fix(Fs / Fbp(1));
    end
    if N > floor( (size(dat,2) - 1) / 3)
      N=floor(size(dat,2)/3) - 1;
    end
    [B, A] = fir1(N, [min(Fbp)/Fn max(Fbp)/Fn]);
  case 'firls' % from NUTMEG's implementation
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
    z = zeros(1,length(f));
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
    z(pos1:pos2) = 1;
    A = 1;
    B = firls(N,f,z); % requires Matlab signal processing toolbox
  otherwise
    error('unsupported filter type "%s"', type);
end

% demean the data before filtering
meandat = mean(dat,2);
dat = bsxfun(@minus, dat, meandat);

try
  filt = filter_with_correction(B,A,dat,dir);
catch
  switch instabilityfix
    case 'no'
      rethrow(lasterror);
    case 'reduce'
      warning('backtrace', 'off')
      warning('instability detected - reducing the %dth order filter to an %dth order filter', N, N-1);
      warning('backtrace', 'on')
      filt = ft_preproc_bandpassfilter(dat,Fs,Fbp,N-1,type,dir,instabilityfix);
    case 'split'
      N1 = ceil(N/2);
      N2 = floor(N/2);
      warning('backtrace', 'off')
      warning('instability detected - splitting the %dth order filter in a sequential %dth and a %dth order filter', N, N1, N2);
      warning('backtrace', 'on')
      filt1 = ft_preproc_bandpassfilter(dat  ,Fs,Fbp,N1,type,dir,instabilityfix);
      filt  = ft_preproc_bandpassfilter(filt1,Fs,Fbp,N2,type,dir,instabilityfix);
    otherwise
      error('incorrect specification of instabilityfix');
  end % switch
end

