function [filt] = lowpassfilter(dat,Fs,Flp,N,type,dir)

% LOWPASSFILTER removes high frequency components from EEG/MEG data
% 
% Use as
%   [filt] = lowpassfilter(dat, Fsample, Flp, N, type, dir)
% where
%   dat        data matrix (Nchans X Ntime)
%   Fsample    sampling frequency in Hz
%   Flp        filter frequency
%   N          optional filter order, default is 6 (but) or 25 (fir)
%   type       optional filter type, can be
%                'but' Butterworth IIR filter (default)
%                'fir' FIR filter using MATLAB fir1 function 
%   dir        optional filter direction, can be
%                'onepass'         forward filter only
%                'onepass-reverse' reverse filter only, i.e. backward in time
%                'twopass'         zero-phase forward and reverse filter (default)
%
% Note that a one- or two-pass filter has consequences for the
% strength of the filter, i.e. a two-pass filter with the same filter
% order will attenuate the signal twice as strong.
%
% See also HIGHPASSFILTER, BANDPASSFILTER

% Copyright (c) 2003, Robert Oostenveld
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

% set the default filter order later
if nargin<4
    N = [];
end

% set the default filter type
if nargin<5
  type = 'but';
end

% set the default filter direction
if nargin<6
  dir = 'twopass';
end

% Nyquist frequency
Fn = Fs/2;

% compute filter coefficients
switch type
  case 'but'
    if isempty(N)
      N = 6;
    end
    [B, A] = butter(N, max(Flp)/Fn);
  case 'fir'
    if isempty(N)
      N = 25;
    end
    [B, A] = fir1(N, max(Flp)/Fn);
end  

% apply filter to the data
switch dir
  case 'onepass'
    filt = filter(B, A, dat')';
  case 'onepass-reverse'
    dat  = fliplr(dat);
    filt = filter(B, A, dat')';
    filt = fliplr(filt);
  case 'twopass'
    filt = filtfilt(B, A, dat')';
end

