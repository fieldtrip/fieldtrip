function [interp] = align_ijk2xyz(interp)

% ALIGN_IJK2XYZ flips and permutes the 3D volume data such that the
% voxel indices and the headcorodinates approximately correspond. The
% homogenous transformation matrix is modified accordingly, to ensure that
% the headcoordinates of each individual voxel do not change.

if isfield(interp, 'inside') && isfield(interp, 'outside')
  % reformat the inside field to a full volume so that it can be flipped/permuted as well
  tmp = [];
  tmp(interp.inside)  = 1;
  tmp(interp.outside) = 0;
  interp.inside  = tmp;
  interp.outside = [];
end

% determine the known volume parameters
param = parameterselection('all', interp);

% ensure that each of the volumes is 3D
for i=1:length(param)
  interp = setsubfield(interp, param{i}, reshape(getsubfield(interp, param{i}), interp.dim));
end

% permute the 3D volume so that the indices and head-coordinate axes correspond approximately
[dum, dim1] = max(abs(interp.transform(1,1:3)));
[dum, dim2] = max(abs(interp.transform(2,1:3)));
[dum, dim3] = max(abs(interp.transform(3,1:3)));
dim = [dim1 dim2 dim3];
if length(unique(dim))<3
  error('could not determine the correspondence between volume and headcoordinate axes');
else
  for i=1:length(param)
    interp = setsubfield(interp, param{i}, permute(getsubfield(interp, param{i}), dim));
  end
  interp.transform(:,1:3) = interp.transform(:,dim);
end

% update the dimensions of the volume
interp.dim = size(getsubfield(interp, param{1}));
if isfield(interp, 'xgrid')
  interp.xgrid = 1:interp.dim(1);
  interp.ygrid = 1:interp.dim(2);
  interp.zgrid = 1:interp.dim(3);
end

% subsequently flip the volume along each direction, so that the diagonal of the transformation matrix is positive
flipx = eye(4); flipx(1,1) = -1; flipx(1,4) = interp.dim(1)+1;
flipy = eye(4); flipy(2,2) = -1; flipy(2,4) = interp.dim(2)+1;
flipz = eye(4); flipz(3,3) = -1; flipz(3,4) = interp.dim(3)+1;
if interp.transform(1,1)<0
  for i=1:length(param)
    interp = setsubfield(interp, param{i}, flipdim(getsubfield(interp, param{i}), 1));
  end
  interp.transform = interp.transform * flipx;
end
if interp.transform(2,2)<0
  for i=1:length(param)
    interp = setsubfield(interp, param{i}, flipdim(getsubfield(interp, param{i}), 2));
  end
  interp.transform = interp.transform * flipy;
end
if interp.transform(3,3)<0
  for i=1:length(param)
    interp = setsubfield(interp, param{i}, flipdim(getsubfield(interp, param{i}), 3));
  end
  interp.transform = interp.transform * flipz;
end

if isfield(interp, 'inside') && isfield(interp, 'outside')
  % restore the inside and outside fields in their proper formats
  tmp = interp.inside(:);
  tmp(find(isnan(tmp))) = 0; % nan values also indicate outside the brain
  interp.inside  = find( tmp);
  interp.outside = find(~tmp);
end

