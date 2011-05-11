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

if ~all(size(transform)==4) || any(transform(4,:) ~= [0 0 0 1])
  error('invalid transformation matrix');
end

if ft_senstype(sens, 'eeg')

  % any normal coordinate transformation is in principle fine
  % apply the translation, rotation and possibly scaling to the electrode positions
  sens.pnt = apply(transform, sens.pnt);

elseif ft_senstype(sens, 'meg')

  % only a rigid body transformation (translation+rotation) without rescaling is allowed
  rotation = eye(4);
  rotation(1:3,1:3) = transform(1:3,1:3);

  if abs(det(rotation)-1)>10*eps
    error('only a rigid body transformation without rescaling is allowed for MEG sensors');
  end

  % apply the translation and rotation to the coil positions
  sens.pnt = apply(transform, sens.pnt);
  % the sensor coil orientations should be rotated but not translated
  sens.ori = apply(rotation, sens.ori);

else
  error('unsupported or unrecognized type of sensors');
end

if isfield(sens, 'fid') && isfield(sens.fid, 'pnt')
  % also apply the translation, rotation and scaling to the fiducial points
  sens.fid.pnt = apply(transform, sens.pnt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that applies the homogenous transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new] = apply(transform, old)
old(:,4) = 1;
new = old * transform';
new = new(:,1:3);
