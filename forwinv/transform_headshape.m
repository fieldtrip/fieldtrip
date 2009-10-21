function [shape] = transform_headshape(transform, shape)

% TRANSFORM_HEADSHAPE applies a homogenous coordinate transformation to a
% structure with headshape and fiducial information.
%
% Use as
%   shape = transform_headshape(transform, shape)
%
% See also READ_HEADSHAPE

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: transform_headshape.m,v $
% Revision 1.3  2008/07/21 20:31:00  roboos
% also support gradiometer coil orientation (only rotate)
%
% Revision 1.2  2008/04/15 15:33:58  roboos
% fixed small bug (thanks to Vladimir)
%
% Revision 1.1  2008/04/11 16:14:03  roboos
% first implementation, simple helper function for spm integration and symmetry with vol and sens
%

if any(transform(4,:) ~= [0 0 0 1])
  error('invalid transformation matrix');
end

if isfield(shape, 'pnt') && ~isempty(shape.pnt)
  % this also works if the structure describes electrode or gradiometer positions instead of a headshape
  shape.pnt = apply(transform, shape.pnt);
end

if isfield(shape, 'ori')
  % gradiometer coil orientations should only be rotated and not translated
  rotation = eye(4);
  rotation(1:3,1:3) = transform(1:3,1:3);
  if abs(det(rotation)-1)>10*eps
    error('only a rigid body transformation without rescaling is allowed for MEG sensors');
  end
  % apply the rotation to the coil orientations
  shape.ori = apply(rotation, sens.ori);
end

if isfield(shape, 'fid') && isfield(shape.fid, 'pnt')
  % apply the same transformation on the fiducials
  shape.fid.pnt = apply(transform, shape.fid.pnt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that applies the homogenous transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new] = apply(transform, old)
old(:,4) = 1;
new = old * transform';
new = new(:,1:3);
