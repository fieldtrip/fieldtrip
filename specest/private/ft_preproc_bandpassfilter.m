function [filt] = ft_preproc_bandpassfilter(dat, Fs, Fbp, N, type, dir)

% FT_PREPROC_BANDPASSFILTER applies a band-pass filter to the data and thereby
% removes the spectral components in the data except for the ones in the
% specified frequency band
%
% Use as
%   [filt] = ft_preproc_bandpassfilter(dat, Fsample, Fbp, N, type, dir)
% where
%   dat        data matrix (Nchans X Ntime)
%   Fsample    sampling frequency in Hz
%   Fbp        frequency band, specified as [Fhp Flp]
%   N          optional filter order, default is 4 (but) or 25 (fir)
%   type       optional filter type, can be
%                'but' Butterworth IIR filter (default)
%                'fir' FIR filter using Matlab fir1 function
%   dir        optional filter direction, can be
%                'onepass'         forward filter only
%                'onepass-reverse' reverse filter only, i.e. backward in time
%                'twopass'         zero-phase forward and reverse filter (default)
%
% Note that a one- or two-pass filter has consequences for the
% strength of the filter, i.e. a two-pass filter with the same filter
% order will attenuate the signal twice as strong.
%
% See also PREPROC

% Copyright (c) 2003-2008, Robert Oostenveld
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
% $Id: ft_preproc_bandpassfilter.m 2455 2010-12-16 15:57:00Z stekla $

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
      N = 25;
    end
    [B, A] = fir1(N, [min(Fbp)/Fn max(Fbp)/Fn]);
end

filt = filter_with_correction(B,A,dat,dir);