function filt = filter_with_correction(B,A,dat,dir,usefftfilt)

% FILTER_WITH_CORRECTION applies a to the data and corrects
% edge-artifacts for one-pass filtering.
%
% Use as
%   [filt] = filter_with_correction(B,A,dat,dir);
% where
%   B,A        filter coefficients
%   dat        data matrix (Nchans X Ntime)
%   dir        optional filter direction, can be
%                'onepass'                   forward filter only
%                'onepass-reverse'           reverse filter only, i.e. backward in time
%                'twopass'                   zero-phase forward and reverse filter (default)
%                'twopass-reverse'           zero-phase reverse and forward filter
%                'twopass-average'           average of the twopass and the twopass-reverse
%                'onepass-zerophase'         zero-phase forward filter with delay compensation (default for firws, linear-phase symmetric FIR only)
%                'onepass-reverse-zerophase' zero-phase reverse filter with delay compensation
%                'onepass-minphase'          minimum-phase converted forward filter (non-linear!, firws only)
%
% Note that a one- or two-pass filter has consequences for the
% strength of the filter, i.e. a two-pass filter with the same filter
% order will attenuate the signal twice as strong.

% Copyright (c) 2010, Stefan Klanke
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

% convert the data to double precision
% see  http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2653
inputclass = class(dat);
B = double(B);
A = double(A);
dat = double(dat);

poles = roots(A);
if any(abs(poles) >= 1)
  ft_error('Calculated filter coefficients have poles on or outside the unit circle and will not be stable. Try a higher cutoff frequency or a different type/order of filter.');
end

dcGain = sum(B)/sum(A);

[d,N] = size(dat);

switch dir
  case 'onepass'
    offset = dat(:,1);
    dat = dat - repmat(offset,1,N);
    filt = filter(B, A, dat')' + repmat(dcGain*offset, 1, N);
  case 'onepass-reverse'
    offset = dat(:,end);
    dat  = fliplr(dat) - repmat(offset,1,N);
    filt = filter(B, A, dat')';
    filt = fliplr(filt) + repmat(dcGain*offset, 1, N);
  case 'twopass'
    % filtfilt does the correction for us
    filt = filtfilt(B, A, dat')';
  case 'twopass-reverse'
    % filtfilt does the correction for us
    filt = fliplr(filtfilt(B, A, fliplr(dat)')');
  case 'twopass-average'
    % take the average from the twopass and the twopass-reverse
    filt1 = filtfilt(B, A, dat')';
    filt2 = fliplr(filtfilt(B, A, fliplr(dat)')');
    filt  = (filt1 + filt2)/2;
  case 'onepass-zerophase'
    filt = fir_filterdcpadded(B, A, dat', 0, usefftfilt)';
  case 'onepass-reverse-zerophase'
    offset = dat(:,end);
    dat  = fliplr(dat) - repmat(offset,1,N);
    filt = fir_filterdcpadded(B, A, dat', 0, usefftfilt)';
    filt = fliplr(filt) + repmat(dcGain*offset, 1, N);
  case 'onepass-minphase'
    filt = fir_filterdcpadded(B, A, dat', 1, usefftfilt)';
  otherwise
    ft_error('unsupported filter direction "%s"', dir);
end

% cast it back into the type of the input data, which can e.g. be single or int32
filt = cast(filt, inputclass);
