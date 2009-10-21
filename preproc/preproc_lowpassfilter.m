function [filt] = preproc_lowpassfilter(dat,Fs,Flp,N,type,dir)

% PREPROC_LOWPASSFILTER applies a low-pass filter to the data and thereby
% removes all high frequency components in the data
% 
% Use as
%   [filt] = preproc_lowpassfilter(dat, Fsample, Flp, N, type, dir)
% where
%   dat        data matrix (Nchans X Ntime)
%   Fsample    sampling frequency in Hz
%   Flp        filter frequency
%   N          optional filter order, default is 6 (but) or 25 (fir)
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
% $Log: preproc_lowpassfilter.m,v $
% Revision 1.2  2008/05/23 09:13:58  roboos
% cleaned up code and documentation, ensure that all functions are consistent, added proper implementation to the scratch functions
%
% Revision 1.1  2008/05/23 06:54:22  roboos
% created initial scratch version of preprocessing module, to be used in fieldtrip or as stand-alone toolbox (e.g. in spm8 or braingain)
% some functions are copies of existing roboos/misc versions, some just contain some example code for the implementation
%
% Revision 1.5  2006/08/31 07:57:22  roboos
% implemented onepass-reverse filter: usefull for stimulus artifacts, e.g. TMS or electrical stimulation
%
% Revision 1.4  2006/06/14 12:36:21  roboos
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
% Revision 1.1  2003/04/04 09:53:37  roberto
% new implementation, using 6th order Butterworth FIR filter
%

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

