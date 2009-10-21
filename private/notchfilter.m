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
% $Log: notchfilter.m,v $
% Revision 1.6  2007/04/16 16:08:52  roboos
% removed the old code that was commented out already
%
% Revision 1.5  2004/01/21 13:04:21  roberto
% changed help and comments
%
% Revision 1.4  2003/06/12 08:40:44  roberto
% added variable option to determine filter order
% changed default order from 6 to 4 for notch and bandpass
%
% Revision 1.3  2003/04/04 09:53:06  roberto
% played around and tested 3 options, finally implemented 6th
% order Butterworth FIR filter
%
% Revision 1.2  2003/03/28 15:07:57  roberto
% added default for Fl (50 Hz)
% fixed bug with channel iterator in for-loop (i, was also used for complex number)
%
% Revision 1.1  2003/03/27 09:03:40  roberto
% new implementation based on input from pascal and Ole
%

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

