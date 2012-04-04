function dat = ft_preproc_medianfilter(dat, order);

% FT_PREPROC_MEDIANFILTER applies a median filter, which smooths the data with
% a boxcar-like kernel except that it keeps steps in the data. This
% function requires the Matlab Signal Processing toolbox.
%
% Use as
%   [dat] = ft_preproc_medianfilter(dat, order)
% where
%   dat        data matrix (Nchans X Ntime)
%   order      number, the length of the median filter kernel (default = 25)
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

% set the default filter order
if nargin<2 || isempty(order)
  error('the order of the median filter is not specified');
end

dat = medfilt1(dat, order, [], 2);
