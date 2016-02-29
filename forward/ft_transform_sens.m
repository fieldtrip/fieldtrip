function [sens] = ft_transform_sens(transform, sens)

% FT_TRANSFORM_SENS applies a homogenous coordinate transformation to a
% structure with EEG electrodes or MEG gradiometers. For MEG gradiometers
% the homogenous transformation matrix should be limited to a rigid-body
% translation plus rotation.
%
% Use as
%   sens = ft_transform_sens(transform, sens)
%
% See also FT_READ_SENS, FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2008-2011, Robert Oostenveld
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

sens = ft_transform_geometry(transform, sens);
