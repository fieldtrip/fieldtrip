function [vol] = ft_transform_vol(transform, vol)

% FT_TRANSFORM_VOL applies a homogenous coordinate transformation to
% a structure with an EEG or MEG colume conduction model. The homogenous
% transformation matrix should be limited to a rigid-body translation
% plus rotation and a global rescaling.
%
% Use as
%   vol = ft_transform_vol(transform, vol)
%
% See also FT_READ_VOL, FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if any(transform(4,:) ~= [0 0 0 1])
  error('invalid transformation matrix');
end

% only a rigid body transformation without rescaling is allowed
rotation = eye(4);
rotation(1:3,1:3) = transform(1:3,1:3);

% FIXME insert check for nonuniform scaling, should give an error

% if abs(det(rotation)-1)>10*eps
%   error('only a rigid body transformation without rescaling is allowed');
% end

switch ft_voltype(vol)
  
  case {'singlesphere' 'multisphere' 'concentric'}
    if isfield(vol, 'o')
      % shift the center of the spheres, an optional rotation does not affect them
      vol.o = apply(transform, vol.o);
    end

  case {'bem', 'dipoli', 'bemcp', 'asa', 'avo', 'nolte'}
    for i=1:length(vol.bnd)
      % apply the transformation to each of the triangulated surface descriptions
      vol.bnd(i).pnt = apply(transform, vol.bnd(i).pnt);
      if isfield(vol.bnd(i), 'nrm')
        % also apply it to the surface normals
        vol.bnd(i).nrm = apply(transform, vol.bnd(i).nrm);
      end
    end

  case 'infinite'
    % nothing to do, since it is an infinite vacuum

  case 'neuromag'
    error('not supported for neuromag forward model');

  otherwise
    error('unsupported or unrecognized type of volume conductor model');
end % switch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that applies the homogenous transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new] = apply(transform, old)
old(:,4) = 1;
new = old * transform';
new = new(:,1:3);
