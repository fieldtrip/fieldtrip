function [vol] = transform_vol(transform, vol)

% TRANSFORM_VOL applies a homogenous coordinate transformation to
% a structure with an EEG or MEG colume conduction model. The homogenous
% transformation matrix should be limited to a rigid-body translation
% plus rotation and a global rescaling.
%
% Use as
%   vol = transform_vol(transform, vol)
%
% See also READ_VOL, PREPARE_VOL_SENS, COMPUTE_LEADFIELD

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: transform_vol.m,v $
% Revision 1.6  2009/02/06 08:31:19  roboos
% added bemcp as volume type
%
% Revision 1.5  2008/04/18 13:16:25  roboos
% removed check for scaling
%
% Revision 1.4  2008/04/15 20:36:21  roboos
% added explicit handling of various BEM implementations, i.e. for all voltype variants
%
% Revision 1.3  2008/03/06 09:27:31  roboos
% updated documentation
%
% Revision 1.2  2008/03/05 15:18:41  roboos
% the previous version was still empty, I now made a proper implementation for the translation of various objects
%
% Revision 1.1  2008/01/28 20:31:28  roboos
% initial implementation, empty stubs
%

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

switch voltype(vol)
  
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
