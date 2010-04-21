function [dat, beta, x] = ft_preproc_detrend(dat, begsample, endsample, order)

% FT_PREPROC_DETREND removes linear or higher order polynomial trends from the
% data using using General Linear Modeling
%
% Use as
%   [dat] = ft_preproc_detrend(dat, begin, end, order)
% where
%   dat        data matrix (Nchans X Ntime)
%   begsample  index of the begin sample for the trend estimate
%   endsample  index of the end sample for the trend estimate
%   order      number representing the polynomial order (default = 1, i.e. linear)
%
% If no begin and end sample are specified for the trend estimate, it
% will be estimated on the complete data.
%
% See also PREPROC

% Copyright (C) 2008, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% determine the size of the data
[Nchans, Nsamples] = size(dat);

% determine the interval to use for baseline correction
if nargin<2 || isempty(begsample)
  begsample = 1;
end
if nargin<3 || isempty(endsample)
  endsample = Nsamples;
end

% determine the order of the polynomial trend to be removed, default is linear
if nargin<4 || isempty(order)
  order = 1;
end

% create a matrix with regressor components
basis    = 1:Nsamples;
x        = zeros(order+1,Nsamples);
for i=0:order
  x(i+1,:) = basis.^(i);
end
% estimate the polynomial trend using General Linear Modeling, where dat=beta*x+noise
beta = dat(:,begsample:endsample)/x(:,begsample:endsample);
% subtract the trend from the complete data
dat  = dat - beta*x;
