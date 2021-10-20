function y = ft_preproc_slidingrange(dat, width, normalize, varargin)

% FT_PREPROC_SLIDINGRANGE computes the range of the data in a sliding time
% window of the width specified.
%
% Use as
%   [dat] = ft_preproc_slidingrange(dat, width, normalize)
% where
%   dat        data matrix (Nchans x Ntime)
%   width      width of the smoothing kernel, this should be an odd number since the window needs to be centered on an individual sample
%   normalize  boolean, whether to normalize the range of the data with the square root of the window size (default = false)
%
% If the data contains NaNs, these are ignored for the computation, but retained in
% the output.
%
% See also PREPROC

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

if nargin>3
  % for backward compatibility, this function was first implemented as ft_preproc_slidingrange(dat, width, ...)
  % with 'normalize' as an optional key-value pair, but all other PREPROC functions take fixed input arguments
  normalize = varargin{1};
end

if nargin<3 || isempty(normalize)
  normalize = false;
end

% preprocessing fails on channels that contain NaN
if any(isnan(dat(:)))
  ft_warning('FieldTrip:dataContainsNaN', 'data contains NaN values');
end

if mod(width+1, 2)
  ft_error('width should be an odd number');
end

% compute half width
h = (width-1)/2;

n = size(dat,2);
minval = zeros(size(dat));
maxval = zeros(size(dat));

for i=1:n
  begsample = i-h;
  endsample = i+h;
  if begsample<1
    begsample = 1;
  end
  if endsample>n
    endsample=n;
  end
  minval(:,i) = min(dat(:,begsample:endsample),[],2);
  maxval(:,i) = max(dat(:,begsample:endsample),[],2);
end

y = maxval - minval;

if istrue(normalize)
  y = y ./ sqrt(width);
end
