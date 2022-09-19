function [R] = ft_connectivity_cancorr(inputdata, varargin)

% FT_CONNECTIVITY_CANCORR computes the canonical correlation or canonical coherence
% between multiple variables. Canonical correlation analysis (CCA) is a way of
% measuring the linear relationship between two multidimensional variables. It finds
% two bases, one for each variable, that are optimal with respect to correlations
% and, at the same time, it finds the corresponding correlations.
%
% Use as
%   [R] = ft_connectivity_cancorr(inputdata, ...)
%
% The input data should be a covariance or cross-spectral density array organized as
%   Channel x Channel
% or
%   Channel x Channel (x Frequency)
%
% The output R represents the max(indices)*max(indices) canonical correlation matrix
% or canonical coherence matrix.
%
% Additional optional input arguments come as key-value pairs:
%   'indices'  = 1xNchan vector with indices of the groups to which the channels belong,
%                e.g. [1 1 2 2] for a 2-by-2 connectivity between 2 planar MEG channels
%   'realflag' = boolean flag whether to use the real-valued part only for the determination
%                of the rotation (default = false)
%
% See also CONNECTIVITY, FT_CONNECTIVITYANALYSIS

% Copyright (C) 2021 Jan-Mathijs Schoffelen
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

indices  = ft_getopt(varargin, 'indices');
realflag = ft_getopt(varargin, 'realflag', false);

if isempty(indices) && isequal(size(inputdata(:,:,1)), [2 2])
  % simply assume two channels
  indices = [1 1 2 2];
end

sizein  = size(inputdata);
sizeout = [sizein 1];
sizeout(1:2) = max(indices);

R = zeros(sizeout);
for kk = 1:sizeout(3)
  for k = 1:sizeout(1)
    for m = 1:sizeout(1)
      indx1 = find(indices==k);
      indx2 = find(indices==m);
      
      [E, D] = multivariate_decomp(inputdata([indx1 indx2], [indx1 indx2], kk), 1:numel(indx1), numel(indx1)+(1:numel(indx2)), 'cca', realflag);
      
      R(k, m, kk) = D(1);
    end
  end
end
