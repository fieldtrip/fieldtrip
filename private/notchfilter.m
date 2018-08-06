function [filt] = notchfilter(dat,Fs,Fl,N)

% NOTCHFILTER line noise reduction filter for EEG/MEG data
%
% [filt] = notchfilter(dat, Fsample, Fline)
%
% where
%   dat        data matrix (Nchans X Ntime)
%   Fsample    sampling frequency in Hz
%   Fline      line noise frequency (would normally be 50Hz)
%   N          optional filter order, default is 4
%
% if Fline is specified as 50, a band of 48-52 is filtered out
% if Fline is specified as [low high], that band is filtered out

% original      (c) 2003, Pascal Fries
% modifications (c) 2003, Robert Oostenveld
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

if nargin<4
  % set the default filter order
  N = 4;
end

Nchans   = size(dat,1);
Nsamples = size(dat,2);

% use a digital FIR filter
Fn = Fs/2;           % Nyquist frequency
if length(Fl)==1
  % default use a notch-width of 2Hz in both directions
  % otherwise use the specified band
  Fl = [Fl-2 Fl+2];
end
[B, A] = butter(N, [min(Fl)/Fn max(Fl)/Fn], 'stop');
filt = filtfilt(B, A, dat')';

