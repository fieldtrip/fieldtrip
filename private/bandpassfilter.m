function [filt] = bandpassfilter(dat,Fs,Fbp,N,type,dir)

% BANDPASSFILTER filters EEG/MEG data in a specified band
% 
% Use as
%   [filt] = bandpassfilter(dat, Fsample, Fbp, N, type, dir)
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
% See also LOWPASSFILTER, HIGHPASSFILTER

% Copyright (c) 2003, Robert Oostenveld
%
% $Log: bandpassfilter.m,v $
% Revision 1.5  2006/08/31 07:57:22  roboos
% implemented onepass-reverse filter: usefull for stimulus artifacts, e.g. TMS or electrical stimulation
%
% Revision 1.4  2006/06/14 12:36:20  roboos
% added the filter direction as additional option, default is 'twopass', i.e. using filtfilt
%
% Revision 1.3  2004/02/11 08:55:13  roberto
% added optional fir1 filter (default still is butterworth), changed
% layout of code for better support of multiple optional arguments,
% extended documentation
%
% Revision 1.2  2003/06/12 08:40:44  roberto
% added variable option to determine filter order
% changed default order from 6 to 4 for notch and bandpass
%
% Revision 1.1  2003/04/04 09:53:36  roberto
% new implementation, using 6th order Butterworth FIR filter
%

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
      N = 4;
    end
    [B, A] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn]);
  case 'fir'
    if isempty(N)
      N = 25;
    end
    [B, A] = fir1(N, [min(Fbp)/Fn max(Fbp)/Fn]);
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

