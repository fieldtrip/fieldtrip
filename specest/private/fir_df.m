function [df, maxDf] = fir_df(cutoffArray, Fs)

% FIR_DF computes default and maximum possible transition band width from
% FIR filter cutoff frequency(ies)
%
% Use as
%   [df, maxDf] = fir_df(cutoffArray, Fs)
% where
%   cutoffArray filter cutoff frequency(ies)
%   Fs          sampling frequency in Hz
%
% Required filter order/transition band width is estimated with the
% following heuristic: transition band width is 25% of the lower cutoff
% frequency, but not lower than 2 Hz, where possible (for bandpass,
% highpass, and bandstop) and distance from passband edge to critical
% frequency (DC, Nyquist) otherwise. 
%
% See also FIRWS, FIRWSORD, INVFIRWSORD

% Copyright (c) 2014, Andreas Widmann
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin < 2 || isempty(cutoffArray) || isempty(Fs)
    ft_error('Not enough input arguments.')
end

% Constants
TRANSWIDTHRATIO = 0.25;
Fn = Fs / 2;

% Max possible transition band width
cutoffArray = sort(cutoffArray);
maxTBWArray = [cutoffArray * 2 (Fn - cutoffArray) * 2 diff(cutoffArray)];
maxDf = min(maxTBWArray);

% Default filter order heuristic
df = min([max([cutoffArray(1) * TRANSWIDTHRATIO 2]) maxDf]);

end
