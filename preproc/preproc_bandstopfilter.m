function [filt] = preproc_bandstopfilter(dat,Fs,Fbp,N,type,dir)

% PREPROC_BANDSTOPFILTER applies a band-stop filter to the data and thereby
% removes the spectral components in the specified frequency band 
% 
% 
% Use as
%   [filt] = preproc_bandstopfilter(dat, Fsample, Fbp, N, type, dir)
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

% Copyright (c) 2007-2008, Robert Oostenveld
%
% $Log: preproc_bandstopfilter.m,v $
% Revision 1.2  2008/05/23 09:13:58  roboos
% cleaned up code and documentation, ensure that all functions are consistent, added proper implementation to the scratch functions
%
% Revision 1.1  2008/05/23 06:54:21  roboos
% created initial scratch version of preprocessing module, to be used in fieldtrip or as stand-alone toolbox (e.g. in spm8 or braingain)
% some functions are copies of existing roboos/misc versions, some just contain some example code for the implementation
%
% Revision 1.1  2007/09/11 15:30:59  roboos
% new implementation based on bandpassfilter
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
      N = 4;
    end
    [B, A] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn], 'stop');
  case 'fir'
    if isempty(N)
      N = 25;
    end
    [B, A] = fir1(N, [min(Fbp)/Fn max(Fbp)/Fn], 'stop');
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

