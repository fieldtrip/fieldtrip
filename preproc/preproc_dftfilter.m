function [filt] = preproc_dftfilter(dat, Fs, Fl)

% PREPROC_DFTFILTER applies a notch filter to the data to remove the 50Hz
% or 60Hz line noise components. This is done by fitting a sine and cosine
% at the specified frequency to the data and subsequently subtracting the
% estimated components. The longer the data is, the sharper the spectral
% notch will be that is removed from the data.
%
% Use as
%   [filt] = preproc_dftfilter(dat, Fsample, Fline)
% where
%   dat        data matrix (Nchans X Ntime)
%   Fsample    sampling frequency in Hz
%   Fline      line noise frequency
%
% The line frequency should be specified as a single number. 
% If omitted, a European default of 50Hz will be assumed.
%
% Preferaby the data should have a length that is a multiple of the 
% oscillation period of the line noise (i.e. 20ms for 50Hz noise). If the
% data is of different lenght, then only the first N complete periods are
% used to estimate the line noise. The estimate is subtracted from the
% complete data.
%
% See also PREPROC

% original      Copyright (C) 2003, Pascal Fries
% modifications Copyright (C) 2003-2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% determine the size of the data
[Nchans, Nsamples] = size(dat);

% set the default filter frequency
if nargin<3 || isempty(Fl)
  Fl = 50;
end

% determine the largest integer number of line-noise cycles that fits in the data
sel = 1:round(floor(Nsamples * Fl/Fs) * Fs/Fl);

% temporarily remove mean to avoid leakage
mdat = mean(dat(:,sel),2);
dat  = dat - mdat(:,ones(1,Nsamples));

% fit a sin and cos to the signal and subtract them
time  = (0:Nsamples-1)/Fs;
tmp  = exp(j*2*pi*Fl*time);                    % complex sin and cos
% ampl = 2*dat*tmp'/Nsamples;                  % estimated amplitude of complex sin and cos
ampl = 2*dat(:,sel)*tmp(sel)'/length(sel);     % estimated amplitude of complex sin and cos on integer number of cycles
est  = ampl*tmp;                               % estimated signal at this frequency
%filt = dat - est;                              % subtract estimated signal
filt = dat - est + mdat(:,ones(1,Nsamples));
filt = real(filt);
