% rcepsminphase() - Convert FIR filter coefficient to minimum phase
%
% Usage:
%   >> b = minphaserceps(b);
%
% Inputs:
%   b - FIR filter coefficients
%
% Outputs:
%   bMinPhase - minimum phase FIR filter coefficients
%
% Author: Andreas Widmann, University of Leipzig, 2013
%
% References:
%   [1] Smith III, O. J. (2007). Introduction to Digital Filters with Audio
%       Applications. W3K Publishing. Retrieved Nov 11 2013, from
%       https://ccrma.stanford.edu/~jos/fp/Matlab_listing_mps_m.html
%   [2] Vetter, K. (2013, Nov 11). Long FIR filters with low latency.
%       Retrieved Nov 11 2013, from
%       http://www.katjaas.nl/minimumphase/minimumphase.html

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2013 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% $Id$

function [bMinPhase] = minphaserceps(b)

% Line vector
b = b(:)';

n = length(b);
upsamplingFactor = 1e3; % Impulse response upsampling/zero padding to reduce time-aliasing
nFFT = 2^ceil(log2(n * upsamplingFactor)); % Power of 2
clipThresh = 1e-8; % -160 dB

% Spectrum
s = abs(fft(b, nFFT));
s(s < clipThresh) = clipThresh; % Clip spectrum to reduce time-aliasing

% Real cepstrum
c = real(ifft(log(s)));

% Fold
c = [c(1) [c(2:nFFT / 2) 0] + conj(c(nFFT:-1:nFFT / 2 + 1)) zeros(1, nFFT / 2 - 1)];

% Minimum phase
bMinPhase = real(ifft(exp(fft(c))));

% Remove zero-padding
bMinPhase = bMinPhase(1:n);

end
