function [headmodel] = ft_transform_vol(transform, headmodel)

% FT_TRANSFORM_VOL applies a homogenous coordinate transformation to
% a structure with an EEG or MEG colume conduction model. The homogenous
% transformation matrix should be limited to a rigid-body translation
% plus rotation and a global rescaling.
%
% Use as
%   headmodel = ft_transform_vol(transform, headmodel)
%
% See also FT_READ_VOL, FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD, FT_TRANSFORM_GEOMETRY

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

headmodel = ft_transform_geometry(transform, headmodel);
