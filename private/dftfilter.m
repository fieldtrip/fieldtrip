function [filt] = dftfilter(dat,Fs,Fl)

% DFTFILTER line noise reduction filter for EEG/MEG data
% 
% [filt] = dftfilter(dat, Fsample, Fline)
%
% where
%   dat        data matrix (Nchans X Ntime)
%   Fsample    sampling frequency in Hz
%   Fline      line noise frequency
%
% The line frequency should be specified as a single number. 
% If omitted, a European default of 50Hz will be assumed.
%
% Preferaby the data should have a length that is a multiple 
% of the period of oscillation of the line noise (20ms for
% 50Hz noise).
%
% See also NOTCHFILTER, 

% original      Copyright (C) 2003, Pascal Fries
% modifications Copyright (C) 2003, Robert Oostenveld
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

if nargin<3 | isempty(Fl)
  Fl = 50;
end

Nchans   = size(dat,1);
Nsamples = size(dat,2);

% determine the largest integer number of line-noise cycles that fits in the data
sel = 1:round(floor(Nsamples * Fl/Fs) * Fs/Fl);

% fit a sin and cos to the signal and subtract them
time  = (0:Nsamples-1)/Fs;
tmp  = exp(j*2*pi*Fl*time);                    % complex sin and cos
% ampl = 2*dat*tmp'/Nsamples;                  % estimated amplitude of complex sin and cos
ampl = 2*dat(:,sel)*tmp(sel)'/length(sel);     % estimated amplitude of complex sin and cos on integer number of cycles
est  = ampl*tmp;                               % estimated signal at this frequency
filt = dat - est;                              % subtract estimated signal
filt = real(filt);

