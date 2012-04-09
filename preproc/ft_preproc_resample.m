function [datout, tim] = ft_preproc_resample(dat, Fold, Fnew, method)

% FT_PREPROC_RESAMPLE resamples all channels in the data matrix
%
% Use as
%   dat = ft_preproc_resample(dat, Fold, Fnew, method)
% where
%   dat    = matrix with the input data (Nchans X Nsamples)
%   Fold   = scalar, original sampling frequency in Hz
%   Fnew   = scalar, desired sampling frequency in Hz
%   method = string, can be 'resample', 'decimate', 'downsample'
%
% The resample method applies an anti-aliasing (lowpass) FIR filter to
% the data during the resampling process, and compensates for the filter's
% delay. For the other two methods you should apply an anti-aliassing
% filter prior to calling this function.
%
% See also PREPROC, FT_PREPROC_LOWPASSFILTER

% Copyright (C) 2006-2012, Robert Oostenveld
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

[nchans, nsamples] = size(dat);

if nargout==2
  tim = 1:size(dat,nsamples);
  tim = ft_preproc_resample(tim, Fold, Fnew, method);
end

if Fold==Fnew
  return
end

typ = class(dat);

% resample and decimate require double formatted input
if ~strcmp(method, 'downsample') && ~strcmp(typ, 'double')
  dat = cast(dat, 'double');
end

switch method
  case 'resample'
    % the actual implementation resamples along columns
    datout = resample(dat', Fnew, Fold)';
    
  case 'decimate'
    fac         = round(Fold/Fnew);
    % this only works one channel at the time
    nresampled  = ceil(nsamples/fac);
    datout      = zeros(nchans, nresampled);
    for i=1:nchans
      datout(i,:) = decimate(dat(i,:), fac);
    end
    
  case 'downsample'
    fac = Fold/Fnew;
    % the actual implementation resamples along columns
    datout = downsample(dat', fac)';
    
  otherwise
    error('unsupported resampling method');
end

if ~strcmp(method, 'downsample') && ~strcmp(typ, 'double')
  % convert back into the original input format
  datout = cast(datout, typ);
end

