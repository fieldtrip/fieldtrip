% evalRationalPolyAscending() - Evaluate rational polynomial as specified
%                               in BVRF file format for mapping numeric values
%                               to measured quantity.
%
% Usage:
%       >> dataout = evalRationalPolyAscending(nv, Num, Denom, RE)
%
% Inputs:
%   NV      - Numeric values from the data file (scalar, vector).
%             These are the raw values contained in the
%             *.bvrd data stream.
%
%   Num     - Numerator coefficient array for the rational polynomial.
%             Must follow the ascending-power.Format: [N0 N1 ... NK].
%
%   Denom   - Denominator coefficient array for the rational polynomial.
%             If empty or missing, the default denominator is 1 (i.e. [1]).
%             Format: [1 D1 ... DL]
%
%   RE      - ResolutionPerBit for this channel. This is the physical scaling
%             factor that converts numeric values to sensor output:
%                   SO = NV * RE
%
% Output:
%   MQ      - Measured Quantity in physical units (as defined by chan.Unit).
%
% Author: Ramon Martinez-Cancino, Brain Products GmbH, 2025
%
% Copyright (C) 2025 Brain Products GmbH
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
function y = evalRationalPolyAscending(nv, num, den, RE)

nv = double(nv);
RE = double(RE);

% (NV * RE) is the polynomial variable
so = nv .* RE;

% Ensure row vectors
num = num(:).';
if isempty(den)
    den = 1;         % default denominator
else
    den = den(:)';
end

% ---- Numerator: N0 + N1*(so) + N2*(so)^2 + ... ----
num_poly = zeros(size(so));
pow = ones(size(so));      % (so)^0
for k = 1:numel(num)
    num_poly = num_poly + num(k) .* pow;
    pow = pow .* so;       % (so)^k
end

% ---- Denominator: 1 + D1*(so) + D2*(so)^2 + ... ----
den_poly = ones(size(so)); % leading 1
pow = so;                  % (so)^1
for k = 2:numel(den)
    den_poly = den_poly + den(k) .* pow;
    pow = pow .* so;
end

y = num_poly ./ den_poly;
end
