function [datout, tim, Fnew] = ft_preproc_resample(dat, Fold, Fnew, method)

% FT_PREPROC_RESAMPLE resamples all channels in the data matrix
%
% Use as
%   dat = ft_preproc_resample(dat, Fold, Fnew, method)
% where
%   dat    = matrix with the input data (Nchans X Nsamples)
%   Fold   = scalar, original sampling frequency in Hz
%   Fnew   = scalar, desired sampling frequency in Hz
%   method = string, can be 'resample', 'decimate', 'downsample', 'fft'
%
% The resample method applies an anti-aliasing (lowpass) FIR filter to
% the data during the resampling process, and compensates for the filter's
% delay. For the other two methods you should apply an anti-aliassing
% filter prior to calling this function.
%
% See also PREPROC, FT_PREPROC_LOWPASSFILTER

% Copyright (C) 2006-2012, Robert Oostenveld
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

[nchans, nsamples] = size(dat);

if nargout>1
  tim = 1:nsamples;
  tim = ft_preproc_resample(tim, Fold, Fnew, method);
end

if Fold==Fnew
  ft_warning('requested sampling rate after resampling is the same as the sampling rate of the data, no resampling is performed');
  datout = dat;
  return
end

% resample and decimate require double formatted input
if ~strcmp(method, 'downsample')
  typ = class(dat);
  dat = cast(dat, 'double');
end

% preprocessing fails on channels that contain NaN
if any(isnan(dat(:)))
  ft_warning('FieldTrip:dataContainsNaN', 'data contains NaN values');
end

switch method
  case 'resample'
    [fold, fnew] = rat(Fold./Fnew); %account for non-integer fs
    Fnew         = Fold.*(fnew./fold); %get new fs exact
    
    % the actual implementation resamples along columns
    datout = resample(dat', fnew, fold)';
    
  case 'decimate'
    fac         = round(Fold/Fnew);
    % this only works one channel at the time
    nresampled  = ceil(nsamples/fac);
    datout      = zeros(nchans, nresampled);
    for i=1:nchans
      datout(i,:) = decimate(dat(i,:), fac);
    end
    
  case 'downsample'
    fac = round(Fold/Fnew);
    % the actual implementation resamples along columns
    datout = downsample(dat', fac)';
    
  case 'fft'
    % Code written for SPM by Jean Daunizeau
    fac         = Fnew/Fold;
    nresampled  = floor(nsamples*fac);
    fac         = nresampled/nsamples;
    datfft      = fftshift(fft(dat,[],2),2);
    middle = floor(size(datfft,2)./2)+1;
    if fac>1 % upsample
      npad = floor((nresampled-nsamples)./2);
      
      if nsamples/2 == floor(nsamples/2)
        datfft(:,1) = []; % throw away non symmetric DFT coef
      end
      
      datfft  = [zeros(size(datfft,1),npad), datfft,zeros(size(datfft,1),npad)];
    else % downsample
      ncut    = floor(nresampled./2);
      datfft  = datfft(:,middle-ncut:middle+ncut);
    end
    datout      = fac*ifft(ifftshift(datfft,2),[],2);
  otherwise
    ft_error('unsupported resampling method');
end

if ~strcmp(method, 'downsample')
  % convert back into the original input format
  datout = cast(datout, typ);
end

