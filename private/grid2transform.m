function [volume] = grid2transform(volume)

% GRID2TRANSFORM ensures that the volume contains a homogenous transformation
% matrix. If needed, a homogenous matrix is constructed from the xgrid/ygrid/zgrid
% fields and those fields are changed to 1:Nx, 1:Ny and 1:Nz
%
% See also TRANSFORM2GRID

if ~isfield(volume, 'xgrid'),     volume.xgrid=1:volume.dim(1); end
if ~isfield(volume, 'ygrid'),     volume.ygrid=1:volume.dim(2); end
if ~isfield(volume, 'zgrid'),     volume.zgrid=1:volume.dim(3); end
if ~isfield(volume, 'transform'), volume.transform=eye(4);      end

% make a local copy for convenience
xgrid = volume.xgrid;
ygrid = volume.ygrid;
zgrid = volume.zgrid;
oldtransform = volume.transform;

% check the input arguments
if any(size(oldtransform)~=[4 4])
  error('incorrect size of homogenous transformation matrix');
elseif ~all(abs(diff(xgrid,2))<=100*eps*max(abs(xgrid)))
  error('xgrid is not monotonously decreasing or increasing');
elseif ~all(abs(diff(ygrid,2))<=100*eps*max(abs(ygrid)))
  error('ygrid is not monotonously decreasing or increasing');
elseif ~all(abs(diff(zgrid,2))<=100*eps*max(abs(zgrid)))
  error('zgrid is not monotonously decreasing or increasing');
end

Nx = length(xgrid);
Ny = length(ygrid);
Nz = length(zgrid);

if any(volume.dim(1:3)~=[Nx Ny Nz])
  error('dimensions do not correspond');
end

xmin = min(xgrid);
ymin = min(ygrid);
zmin = min(zgrid);

xmax = max(xgrid);
ymax = max(ygrid);
zmax = max(zgrid);

if all([xmin ymin zmin]==[1 1 1]) && all([xmax ymax zmax]==volume.dim(1:3))
  % nothing needs to be done, since the homogenous transformation matrix already
  % corresponds to a conversion from voxel indices to world coordinates
else
  fprintf('updating homogenous coordinate transformation matrix\n');
  % three transformations are needed to match the old and new transformation matrices
  % oldtransform = newtransform * t3 * t2 * t1

  % first shift to [0 0 0]
  t1 = eye(4);
  t1(1,4) = -xmin;
  t1(2,4) = -ymin;
  t1(3,4) = -zmin;
  % then scale to a length equal to the new number of voxels minus one
  t2 = eye(4);
  t2(1,1) = (Nx-1)/(xmax-xmin);
  t2(2,2) = (Ny-1)/(ymax-ymin);
  t2(3,3) = (Nz-1)/(zmax-zmin);
  % and shift back to [1 1 1]
  t3 = eye(4);
  t3(1,4) = 1;
  t3(2,4) = 1;
  t3(3,4) = 1;

  % combine all the homogenous transformations in such a way that the new
  % transformation behaves identically on the new x/y/zgrid as the old
  % transformation did on the old x/y/zgrid
  newtransform = oldtransform / (t3 * t2 * t1);

  % update the definition of the cardinal axes of the volume
  % volume.xgrid = 1:Nx;
  % volume.ygrid = 1:Ny;
  % volume.zgrid = 1:Nz;

  % update the homogenous transformation matrix
  volume.transform = newtransform;
end

% remove the cardinal axes of the volume, since they have become trivial
volume = rmfield(volume, 'xgrid');
volume = rmfield(volume, 'ygrid');
volume = rmfield(volume, 'zgrid');

