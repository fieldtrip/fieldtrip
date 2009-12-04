function [sens] = transform_sens(transform, sens)

% TRANSFORM_SENS applies a homogenous coordinate transformation to a
% structure with EEG electrodes or MEG gradiometers. For MEG gradiometers
% the homogenous transformation matrix should be limited to a rigid-body
% translation plus rotation.
%
% Use as
%   sens = transform_sens(transform, sens)
%
% See also READ_SENS, PREPARE_VOL_SENS, COMPUTE_LEADFIELD

% Copyright (C) 2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if any(transform(4,:) ~= [0 0 0 1])
  error('invalid transformation matrix');
end

if senstype(sens, 'eeg')

  % any normal coordinate transformation is in principle fine
  % apply the translation, rotation and possibly scaling to the electrode positions
  sens.pnt = apply(transform, sens.pnt);

elseif senstype(sens, 'meg')

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that applies the homogenous transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new] = apply(transform, old)
old(:,4) = 1;
new = old * transform';
new = new(:,1:3);
