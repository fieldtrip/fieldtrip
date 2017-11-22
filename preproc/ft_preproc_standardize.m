function [dat, state] = ft_preproc_standardize(dat, begsample, endsample, state)

% FT_PREPROC_STANDARDIZE performs a z-transformation or standardization
% of the data. The standardized data will have a zero-mean and a unit
% standard deviation.
%
% Use as
%   [dat] = ft_preproc_standardize(dat, begsample, endsample)
% where
%   dat        data matrix (Nchans dat Ntime)
%   begsample  index of the begin sample for the mean and stdev estimate
%   endsample  index of the end sample for the mean and stdev estimate
%
% If no begin and end sample are specified, it will be estimated on the
% complete data.
%
% See also PREPROC

% Copyright (C) 2008, Robert Oostenveld
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

if nargin<2 || isempty(begsample)
  begsample = 1;
end

if nargin<3 || isempty(endsample)
  endsample = size(dat,2);
end

if nargin<4
  state = [];
end

% preprocessing fails on channels that contain NaN
if any(isnan(dat(:)))
  ft_warning('FieldTrip:dataContainsNaN', 'data contains NaN values');
end

% get the data selection
y = dat(:,begsample:endsample);

% determine the size of the selected data: nChans dat nSamples
[m, n] = size(y);

% compute the sum and sum of squares
s  = sum(y,2);
ss = sum(y.^2,2);

% include the state information from the previous calls
if ~isempty(state)
  s  = s  + state.s;
  ss = ss + state.ss;
  n  = n  + state.n;
end

% compute the mean and standard deviation
my = s ./ n;
sy = sqrt((ss - (s.^2)./n) ./ (n-1));

% standardize the complete input data
dat = (dat - repmat(my, 1, size(dat, 2))) ./ repmat(sy, 1, size(dat, 2));

% remember the state
state.s  = s;  % sum
state.ss = ss; % sum of sqares
state.n  = n;  % number of samples
