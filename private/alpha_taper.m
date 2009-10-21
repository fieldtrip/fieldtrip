function [tap] = alpha_taper(n, f);

% ALPHA_TAPER returns an asymmetric taper that can be used to construct a
% complex wavelet with the peak at a distance of 0.8 times the cycle length
% from the end.
%
% Use as
%   tap = alpha_taper(n, f)
% where
%   n = number of samples
%   f = frequency of desired wavelet, relative to the sampling frequency
%
% The taper will be sufficiently long for a wavelet when n>=5/f.
%
% Example:
%   f = 0.01; % 10 Hz wavelet at 1000 Hz sampling rate
%   plot(alpha_taper(5/f, f)); hold on
%   plot(alpha_taper(5/f, f) .* cos(2*pi*10*(-499:0)/1000), 'r');
%   plot(alpha_taper(5/f, f) .* sin(2*pi*10*(-499:0)/1000), 'g');
%
% This function implements equation 3 from Mitchell, Baker and Baker (2007);
% Muscle Responses to Transcranial Stimulation Depend on Background Oscillatory
% Activity. http://jp.physoc.org/cgi/content/abstract/jphysiol.2007.134031v1
%
% The original paper contains a typo. The equation 3 in the paper reads
%   W(F,t) = -(5/4)*F*t * exp( (1+(5/4)*F*t) * i*2*pi*F*t )
% but should read
%   W(F,t) = -(5/4)*F*t * exp( (1+(5/4)*F*t) + i*2*pi*F*t )
% since then it is equal to
%   W(F,t) = -(5/4)*F*t * exp(1+(5/4)*F*t) * exp(i*2*pi*F*t)
% which is simply
%   W(F,t) = taper(F,t) * exp(i*2*pi*F*t)

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: alpha_taper.m,v $
% Revision 1.2  2007/08/06 15:00:03  roboos
% transposed the output, updated documentation
%
% Revision 1.1  2007/08/01 16:20:27  roboos
% first implementation
%

% time axis expressed in cycles of the desired wavelet frequency
t   = ((-n+1):0) * f;
tap = -(5/4).*t .* exp(1+(5/4).*t);
tap = tap';
