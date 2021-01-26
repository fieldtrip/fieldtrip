function [dat] = ft_preproc_rectify(dat)

% FT_PREPROC_RECTIFY rectifies the data, i.e. converts all samples with a
% negative value into the similar magnitude positive value
%
% Use as
%   [dat] = ft_preproc_rectify(dat)
% where
%   dat        data matrix (Nchans X Ntime)
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

dat = abs(dat);

