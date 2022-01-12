% Copyright (C) 1999 Paul Kienzle
% Copyright (C) 2003 Doug Stewart <dastew@sympatico.ca>
% Copyright (C) 2011 Alexander Klein <alexander.klein@math.uni-giessen.de>
% Copyright (C) 2018 John W. Eaton
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.

% Generate a butterworth filter.
% Default is a discrete space (Z) filter.
%
% [b,a] = butter(n, Wc)
%    low pass filter with cutoff pi*Wc radians
%
% [b,a] = butter(n, Wc, 'high')
%    high pass filter with cutoff pi*Wc radians
%
% [b,a] = butter(n, [Wl, Wh])
%    band pass filter with edges pi*Wl and pi*Wh radians
%
% [b,a] = butter(n, [Wl, Wh], 'stop')
%    band reject filter with edges pi*Wl and pi*Wh radians
%
% [z,p,g] = butter(...)
%    return filter as zero-pole-gain rather than coefficients of the
%    numerator and denominator polynomials.
%
% [...] = butter(...,'s')
%     return a Laplace space filter, W can be larger than 1.
%
% [a,b,c,d] = butter(...)
%  return  state-space matrices
%
% References:
%
% Proakis & Manolakis (1992). Digital Signal Processing. New York:
% Macmillan Publishing Company.

% Author: Paul Kienzle <pkienzle@user.sf.net>
% Modified by: Doug Stewart <dastew@sympatico.ca> Feb, 2003

function [a, b, c, d] = butter(n, W, varargin)

if (nargin>4 || nargin<2) || (nargout>4 || nargout<2)
  error('usage: [b, a] or [z, p, g] or [a,b,c,d] = butter (n, W [, "ftype"][,"s"])');
end

% interpret the input parameters
if (~(length(n)==1 && n == round(n) && n > 0))
  error('butter: filter order n must be a positive integer');
end

stop    = 0;
digital = 1;
for i=1:length(varargin)
  switch varargin{i}
    case 's', digital = 0;
    case 'z', digital = 1;
    case { 'high', 'stop' }, stop = 1;
    case { 'low',  'pass', 'bandpass' }, stop = 0;
    otherwise,  error ('butter: expected [high|stop] or [s|z]');
  end
end

[r, c] = size(W);
if ~(length(W)<=2 && (r==1 || c==1))
  error('butter: frequency must be given as w0 or [w0, w1]');
elseif ~(length(W)==1 || length(W) == 2)
  error ('butter: only one filter band allowed');
elseif length(W)==2 && ~(W(1) < W(2))
  error('butter: first band edge must be smaller than second');
end

if digital && ~all(W >= 0 & W <= 1)
  error('butter: critical frequencies must be in (0 1)');
elseif ~digital && ~all(W >= 0)
  error('butter: critical frequencies must be in (0 inf)');
end

% Prewarp to the band edges to s plane
if digital
  T = 2;       % sampling frequency of 2 Hz
  W = 2/T*tan(pi*W/T);
end

% Generate splane poles for the prototype Butterworth filter
% source: Kuc
C = 1; % default cutoff frequency
pole = C*exp(1i*pi*(2*(1:n) + n - 1)/(2*n));
if mod(n,2) == 1 
  pole((n+1)/2) = -1; 
end  % pure real value at exp(i*pi)
zero = [];
gain = C^n;

% splane frequency transform
[zero, pole, gain] = sftrans(zero, pole, gain, W, stop);

% Use bilinear transform to convert poles to the z plane
if digital
  [zero, pole, gain] = bilinear(zero, pole, gain, T);
end
  
% convert to the correct output form
if nargout<=2
  a = real(gain*poly(zero));
  b = real(poly(pole));
elseif nargout==3
  a = zero(:);
  b = pole(:);
  c = gain;
else
  % output ss results
  [a, b, c, d] = zp2ss (zero, pole, gain);
end
