function [dat] = ft_preproc_derivative(dat, order)

% FT_PREPROC_DERIVATIVE computes the temporal Nth order derivative of the
% data
%
% Use as
%   [dat] = ft_preproc_derivative(dat, order)
% where
%   dat        data matrix (Nchans X Ntime)
%   order      number representing the Nth derivative (default = 1)
%
% If the data contains NaNs, these are ignored for the computation, but
% retained in the output.
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

% set the defaults if options are not specified
if nargin<2 || isempty(order)
  order = 1;
end

% preprocessing fails on channels that contain NaN
if any(isnan(dat(:)))
  ft_warning('FieldTrip:dataContainsNaN', 'data contains NaN values');
end

% compute the derivative
for i=1:order
  dat = gradient(dat);
end
