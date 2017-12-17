%firws() - Designs windowed sinc type I linear phase FIR filter
%
% Usage:
%   >> b = firws(m, f);
%   >> b = firws(m, f, w);
%   >> b = firws(m, f, t);
%   >> b = firws(m, f, t, w);
%
% Inputs:
%   m - filter order (mandatory even)
%   f - vector or scalar of cutoff frequency/ies (-6 dB;
%       pi rad / sample)
%
% Optional inputs:
%   w - vector of length m + 1 defining window {default blackman}
%   t - 'high' for highpass, 'stop' for bandstop filter {default low-/
%       bandpass}
%
% Output:
%   b - filter coefficients
%
% Example:
%   fs = 500; cutoff = 0.5; df = 1;
%   m  = firwsord('hamming', fs, df);
%   b  = firws(m, cutoff / (fs / 2), 'high', windows('hamming', m + 1)); 
%
% References:
%   Smith, S. W. (1999). The scientist and engineer's guide to digital
%   signal processing (2nd ed.). San Diego, CA: California Technical
%   Publishing.
%
% Author: Andreas Widmann, University of Leipzig, 2005
%
% See also:
%   firwsord, invfirwsord, kaiserbeta, windows

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

function [b, a] = firws(m, f, t, w)

    a = 1;

    if nargin < 2
        ft_error('Not enough input arguments');
    end
    if length(m) > 1 || ~isnumeric(m) || ~isreal(m) || mod(m, 2) ~= 0 || m < 2
        ft_error('Filter order must be a real, even, positive integer.');
    end
    f = f / 2;
    if any(f <= 0) || any(f >= 0.5)
        ft_error('Frequencies must fall in range between 0 and 1.');
    end
    if nargin < 3 || isempty(t)
        t = '';
    end
    if nargin < 4 || isempty(w)
        if ~isempty(t) && ~ischar(t)
            w = t;
            t = '';
        else
            w = windows('blackman', (m + 1));
        end
    end
    w = w(:)'; % Make window row vector

    b = fkernel(m, f(1), w);

    if length(f) == 1 && strcmpi(t, 'high')
        b = fspecinv(b);
    end

    if length(f) == 2
        b = b + fspecinv(fkernel(m, f(2), w));
        if isempty(t) || ~strcmpi(t, 'stop')
            b = fspecinv(b);
        end
    end

% Compute filter kernel
function b = fkernel(m, f, w)
    m = -m / 2 : m / 2;
    b(m == 0) = 2 * pi * f; % No division by zero
    b(m ~= 0) = sin(2 * pi * f * m(m ~= 0)) ./ m(m ~= 0); % Sinc
    b = b .* w; % Window
    b = b / sum(b); % Normalization to unity gain at DC

% Spectral inversion
function b = fspecinv(b)
    b = -b;
    b(1, (length(b) - 1) / 2 + 1) = b(1, (length(b) - 1) / 2 + 1) + 1;
