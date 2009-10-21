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
% $Log: dftfilter.m,v $
% Revision 1.6  2005/01/27 17:06:22  roboos
% fixed bug in normalization of sine and cosine amplitude estimate in case number of samples in the data does not match with an integer number of cycles
%
% Revision 1.5  2004/11/17 09:00:01  roboos
% added selection of data to ensure that the sine wave is estimated on an integer number of line-noise cycles
% all data is filtered, only amplitude estimation is done on this selection
%
% Revision 1.4  2003/12/01 08:47:59  roberto
% updated copyright statement
%
% Revision 1.3  2003/10/01 08:46:25  roberto
% updated help
%
% Revision 1.2  2003/10/01 08:45:23  roberto
% updated help
%
% Revision 1.1  2003/10/01 08:45:03  roberto
% first implementation as separate function, used to be notchfilter
%

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

