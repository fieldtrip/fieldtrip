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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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

