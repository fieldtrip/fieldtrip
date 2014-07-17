% plotfresp() - Plot a filter's impulse, step, magnitude, and phase response
%
% Usage:
%   >> plotfresp(b, a, nfft, fs, causal);
%
% Inputs:
%   b     - vector numerator coefficients
%
% Optional inputs:
%   a     - scalar or vector denominator coefficients (IIR support is
%           experimental!) {default 1}
%   nfft  - scalar number of points {default 512}
%   fs    - scalar sampling frequency {default 1}
%   dir   - string filter direction {default 'onepass'}
%
% Author: Andreas Widmann, University of Leipzig, 2005
%
% See also:
%   firws

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005-2014 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

function plotfresp(b, a, nfft, fs, dir)

if nargin < 5 || isempty(dir)
    dir = 'onepass';
end
if nargin < 4  || isempty(fs)
    fs = 1;
end
if nargin < 3 || isempty(nfft)
    nfft = 512;
end
if nargin < 2 || isempty(a)
    a = 1;
end
if nargin < 1
    error('Not enough input arguments.');
end

% FIR?
if isscalar(a) && a == 1
    isFIR = true;
else
    isFIR = false;
end
    
% Linear phase FIR
if isFIR && all(b(:)' == fliplr(b(:)')) % TODO: antisymmetric
    isLinPhaseFir = true;
else
    isLinPhaseFir = false;
end

% Twopass/zerophase?
if strncmp('twopass', dir, 7)
    isTwopass = true;
    isZerophase = true;
elseif strcmp('onepass-zerophase', dir)
    if ~isLinPhaseFir
        error('Onepass-zerophase filtering is only allowed for linear-phase FIR filters.')
    end
    isTwopass = false;
    isZerophase = true;
else
    isTwopass = false;
    isZerophase = false;
end

% Impulse response
if isFIR
    impresp = b(:)';
else
    if ~exist('impz', 'file')
        warning('Plotting IIR filter responses requires signal processing toolbox.')
        return
    end
    impresp = impz(b, a)';
end

% Twopass
if isTwopass
    impresp = conv(impresp, fliplr(impresp));
end
n = length(impresp);

% Zerophase
if isZerophase
    groupdelay = (n - 1) / 2;
    x = -groupdelay:groupdelay;
else
    x = 0:n - 1;
end

nfft = max([2^ceil(log2(n)) nfft]); % Do not truncate impulse response
f = linspace(0, fs / 2, nfft / 2 + 1);
z = fft(impresp, nfft);
z = z(1:nfft / 2 + 1);

% Find open figure window
H = findobj('Tag', 'plotfiltresp', 'type', 'figure');
if ~isempty(H)
    figure(H);
else
    H = figure;
    set(H, 'Tag', 'plotfiltresp');
    posArray = get(H, 'Position');
    posArray(3) = posArray(4) * 1.6;
    set(H, 'Position', posArray);
end

% Formatting
titlePropArray = {'Fontweight', 'bold'};
axisPropArray = {'NextPlot', 'add', 'XGrid', 'on', 'YGrid', 'on', 'Box', 'on'};

% Impulse resonse
ax(1) = subplot(2, 3, 1, axisPropArray{:});
stem(x, impresp, 'fill')
title('Impulse response', titlePropArray{:});
ylabel('Amplitude');

% Step response
ax(4) = subplot(2, 3, 4, axisPropArray{:});
stem(x, cumsum(impresp), 'fill');
title('Step response', titlePropArray{:});
ylimArray = ylim;
if ylimArray(2) < -ylimArray(1) + 1;
    ylimArray(2) = -ylimArray(1) + 1;
    ylim(ylimArray);
end
xMin = []; xMax = [];
childrenArray = get(ax(4), 'Children');
for iChild =1:length(childrenArray)
    xData = get(childrenArray(iChild), 'XData');
    xMin = min([xMin min(xData)]);
    xMax = max([xMax max(xData)]);
end
set(ax([1 4]), 'XLim', [xMin xMax]);
ylabel('Amplitude');

% Magnitude response
ax(2) = subplot(2, 3, 2, axisPropArray{:});
plot(f, abs(z));
title('Magnitude response', titlePropArray{:});
ylabel('Magnitude (linear)');

ax(5) = subplot(2, 3, 5, axisPropArray{:});
plot(f, 20 * log10(abs(z)));
title('Magnitude response', titlePropArray{:});
ylimArray = ylim;
if ylimArray(1) < -200
    ylimArray(1) = -200;
    ylim(ylimArray);
end
ylabel('Magnitude (dB)');

% Phase response
ax(3) = subplot(2, 3, 3, axisPropArray{:});
phaseresp = unwrap(angle(z));
if isZerophase % Correct delay for zero-phase FIR filter?
    delay = -f / fs * groupdelay * 2 * pi;
    phaseresp = phaseresp - delay;
    phaseresp = mod(round(phaseresp / pi), 2) * pi; % Avoid rounding errors; linear-phase FIR only!
end
plot(f, phaseresp);
title('Phase response', titlePropArray{:});
ylabel('Phase (rad)');

% Formatting
xlabelArray = get(ax(1:5), 'XLabel');
if fs == 1
    set([xlabelArray{[2 3 5]}], 'String', 'Normalized frequency (2 \pi rad / sample)');
else
    set([xlabelArray{[2 3 5]}], 'String', 'Frequency (Hz)');
end
set([xlabelArray{[1 4]}], 'String', 'n (samples)');
set(ax([2 3 5]), 'XLim', [0 fs / 2]);
set(ax(1:5), 'ColorOrder', circshift(get(ax(1), 'ColorOrder'), -1));

end
