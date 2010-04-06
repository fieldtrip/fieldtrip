function [dat, tim] = preproc_resample(dat, Fold, Fnew, method)

% PREPROC_RESAMPLE resamples all channels in the data matrix
%
% Use as
%   dat = preproc_resample(dat, Fold, Fnew, method)
% where
%   dat    = matrix with the input data (Nchans X Nsamples)
%   Fold   = scalar, original sampling frequency in Hz
%   Fnew   = scalar, desired sampling frequency in Hz
%   method = string, can be 'resample', 'decimate', 'downsample'

% Copyright (C) 2006-2010, Robert Oostenveld

if nargout==2
  tim = 1:size(dat,2);
  tim = preproc_resample(tim, Fold, Fnew, method);
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

