function [p] = ft_connectivity_plm(input, varargin)

% FT_CONNECTIVITY_PLM computes the phase linearity measurement from a cell
% array of time-domain data, where each cell is an epoch
%
% Use as
%   [p] = ft_connectivity_plm(input, ...)
%
% The input data input should be organized as a cell-array of nchan x ntime signals
%
% Additional optional input arguments come as key-value pairs:
%   bandwidth	=	scalar, half-bandwidth parameter: the frequency range
%			across which to integrate
%   fsample     =       sampling frequency, needed to convert bandwidth to number of bins
%
% The output p contains the phase slope index, v is a variance estimate
% which only can be computed if the data contains leave-one-out samples,
% and n is the number of repetitions in the input data. If the phase slope
% index is positive, then the first chan (1st dim) becomes more lagged (or
% less leading) with higher frequency, indicating that it is causally
% driven by the second channel (2nd dim)
%
% See also FT_CONNECTIVITYANALYSIS

% Copyright (C) 2018, Fabio Baselice, Pierpaolo Sorrentino, Jan-Mathijs Schoffelen
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

% the sequence of steps is as follows:
%  - Hilbert transformation
%  - convolve with complex conjugate
%  - fft
%  - convert bandwidth parameter to number of bins
%  - integrate over bandwidth


% NOTE BY JM: if the user inputs data with different length trials, the fft per trial is going 
% to have different frequency resolutions, which is not good. Better to throw an error in that 
% case.
nsmp = cellfun('size', input, 2);
assert(all(nsmp==nsmp(1)), 'currently there is no support for input, where the trials are of different length'); 

for k = 1:numel(input)
  input{k} = hilbert(input{k}')';
end

% the rest continues here

