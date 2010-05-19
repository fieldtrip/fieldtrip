function [dat, tim] = ft_preproc_resample(dat, Fold, Fnew, method)

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

% Copyright (C) 2006-2010, Robert Oostenveld
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

if nargout==2
  tim = 1:size(dat,2);
  tim = ft_preproc_resample(tim, Fold, Fnew, method);
end

if Fold==Fnew
  return
end

switch method
  case 'resample'
    if ~isa(dat, 'double')
      typ = class(dat);
      dat = typecast(dat, 'double');
      dat = resample(dat, Fnew, Fold);     % this requires a double array
      dat = typecast(dat, typ);
    else
      dat = resample(dat, Fnew, Fold);
    end
  case 'decimate'
    fac = round(Fold/Fnew);
    if ~isa(dat, 'double')
      typ = class(dat);
      dat = typecast(dat, 'double');
      dat = decimate(dat, fac);    % this requires a double array
      dat = typecast(dat, typ);
    else
      dat = decimate(dat, fac);    % this requires a double array
    end
  case 'downsample'
    fac = Fold/Fnew;
    dat = downsample(dat, fac);
  otherwise
    error('unsupported resampling method');
end

