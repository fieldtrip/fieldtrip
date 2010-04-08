function [volume] = transform2grid(volume);

% TRANSFORM2GRID ensures that the volume contains the definition of the
% cardian axes, i.e. xgrid/ygrid/zgrid. If the voluyme contains a
% homogenous coordinate transformation axis that is unequal to eye(4), it
% will try to construct the cardinal axis from that transformation matrix.
%
% See also GRID2TRANSFORM

if ~isfield(volume, 'xgrid'),     volume.xgrid=1:volume.dim(1); end
if ~isfield(volume, 'ygrid'),     volume.ygrid=1:volume.dim(2); end
if ~isfield(volume, 'zgrid'),     volume.zgrid=1:volume.dim(3); end
if ~isfield(volume, 'transform'), volume.transform=eye(4);      end

if all(all(volume.transform==eye(4)))
  % nothing needs to be done, since the homogenous transformation matrix already
  % corresponds to the identity matrix
else
  fprintf('updating values along the cardinal axes\n');
  % check whether it is possible to convert the homogenous transformation
  % matrix into the cardinal axes that are described with xgrid/ygrid/zgrid
  dum = volume.transform .* [
    0 1 1 0
    1 0 1 0
    1 1 0 0
    1 1 1 0 ];
  if any(dum(:)~=0)
    error('rotated coordinate system, cannot compute xgrid/ygrid/zgrid');
  end

  % apply the homogenous coordinate transformation to each of the cardinal axes
  x = (volume.transform * ([volume.xgrid(:) zeros(volume.dim(1),1) zeros(volume.dim(1),1) ones(volume.dim(1),1)]'))';
  x = x(:,1);
  y = (volume.transform * ([zeros(volume.dim(2),1) volume.ygrid(:) zeros(volume.dim(2),1) ones(volume.dim(2),1)]'))';
  y = y(:,2);
  z = (volume.transform * ([zeros(volume.dim(3),1) zeros(volume.dim(3),1) volume.zgrid(:) ones(volume.dim(3),1)]'))';
  z = z(:,3);

  % update the definition of the cardinal axes of the volume
  volume.xgrid = x;
  volume.ygrid = y;
  volume.zgrid = z;

  % and update the homogenous transformation matrix
  % volume.transform = eye(4);
end

% remove the homogenous transformation matrix, since it has become trivial
volume = rmfield(volume, 'transform');

