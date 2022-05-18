% fir_filterdcpadded() - Pad data with DC constant and filter
%
% Usage:
%   >> data = fir_filterdcpadded(b, a, data, causal);
%
% Inputs:
%   b             - vector of filter coefficients
%   a             - 1
%   data          - raw data (times x chans)
%   causal        - boolean perform causal filtering {default 0}
%   usefftfilt    - boolean use fftfilt instead of filter
%
% Outputs:
%   data      - smoothed data
%
% Note:
%   firfiltdcpadded always operates (pads, filters) along first dimension.
%   Not memory optimized.
%
% Author: Andreas Widmann, University of Leipzig, 2014

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

function [ data ] = fir_filterdcpadded(b, a, data, causal, usefftfilt)

% Defaults
if nargin <= 4 || isempty(usefftfilt)
  usefftfilt = 0;
end
if nargin <= 3 || isempty(causal)
  causal = 0;
end

% Check arguments
if nargin < 3
  ft_error('Not enough input arguments.');
end

% Is FIR?
if ~isscalar(a) || a ~= 1
  ft_error('Not a FIR filter. onepass-zerophase and onepass-minphase filtering is available for FIR filters only.')
end

% Group delay
if mod(length(b), 2) ~= 1
  ft_error('Filter order must be a real, even, positive integer.');
end
groupDelay = (length(b) - 1) / 2;

if causal == 0
  % Filter symmetry, add a check with some numeric tolerance, needed if the
  % filter coefficients are computed using the compat functions
  isSym     = all(b(1:groupDelay)  ==  b(end:-1:groupDelay + 2));
  isAntisym = all([b(1:groupDelay) == -b(end:-1:groupDelay + 2) b(groupDelay + 1) == 0]);
  
  % the filter should be either symmetric or antisymmetric: the check above
  % requires exact numeric anti-/symmetry. The compat drop in replacements
  % have some numeric inaccuracies, close to eps, so we might be a bit
  % lenient here
  if ~isSym && isalmostequal(b(1:groupDelay), b(end:-1:groupDelay + 2), 'abstol', eps*10)
    ft_warning('Filter coefficients are nearly symmetric, making them explicitly symmetric');
    bold = b;
    b(1:groupDelay)          = ( bold(1:groupDelay) + bold(end:-1:groupDelay + 2) )./2;
    b(end:-1:groupDelay + 2) = ( bold(1:groupDelay) + bold(end:-1:groupDelay + 2) )./2;
  end
  isSym     = all(b(1:groupDelay)  ==  b(end:-1:groupDelay + 2));
  
  if ~isAntisym && isalmostequal(b(1:groupDelay), -b(end:-1:groupDelay + 2), 'abstol', eps*10) && isalmostequal(b(groupDelay + 1), 0, 'abstol', eps*10)
    ft_warning('Filter coefficients are nearly antisymmetric, making them explicitly antisymmetric');
    bold = b;
    b(1:groupDelay)          =  ( bold(1:groupDelay) - bold(end:-1:groupDelay + 2) )./2;
    b(end:-1:groupDelay + 2) = -( bold(1:groupDelay) - bold(end:-1:groupDelay + 2) )./2;
  end
  isAntisym = all([b(1:groupDelay) == -b(end:-1:groupDelay + 2) b(groupDelay + 1) == 0]);
  
  if ~(isSym || isAntisym)
    ft_error('Filter is not anti-/symmetric. For onepass-zerophase filtering the filter must be anti-/symmetric.')
  end
end

% Padding
if causal
  startPad = repmat(data(1, :), [2 * groupDelay 1]);
  endPad = [];
else
  startPad = repmat(data(1, :), [groupDelay 1]);
  endPad = repmat(data(end, :), [groupDelay 1]);
end

% Filter data (with double precision)
isSingle = isa(data, 'single');

if usefftfilt
  data = fftfilt(double(b), double([startPad; data; endPad]));
else
  data = filter(double(b), 1, double([startPad; data; endPad])); % Pad and filter with double precision
end

% Convert to single
if isSingle
  data = single(data);
end

% Remove padded data
data = data(2 * groupDelay + 1:end, :);

end
