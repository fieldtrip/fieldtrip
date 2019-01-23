function [shape] = ft_transform_headshape(transform, shape)

% FT_TRANSFORM_HEADSHAPE applies a homogenous coordinate transformation to a
% structure with headshape and fiducial information.
%
% Use as
%   shape = ft_transform_headshape(transform, shape)
%
% See also FT_READ_HEADSHAPE, FT_TRANSFORM_GEOMETRY

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

shape = ft_transform_geometry(transform, shape);
