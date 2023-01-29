function [volume, permutevec, flipflags, transform] = align_ijk2xyz(volume)

% ALIGN_IJK2XYZ flips and permutes the 3D volume data such that the axes of
% the voxel indices and the headcoordinates approximately correspond. The
% homogeneous transformation matrix is modified accordingly, to ensure that
% the headcoordinates of each individual voxel do not change. The intention
% is to create a volume structure that has a transform matrix which is
% approximately diagonal in the rotation part.
%
% First, the volume is permuted in order to get the largest (absolute)
% values on the diagonal of the transformation matrix. This permutation is
% reflected by the second output argument.
% Second, the volumes are flipped along the dimensions for which the main
% diagonal elements of the transformation matrix are negative. This is
% reflected by the third output argument.
%
% The second and third argument are in the output in order to be able to
% reverse the operation. Note that in such case first the data have to be
% 'unflipped', and then 'unpermuted' (using ipermute, rather than permute).

% Copyright (C) 2012-2022, Jan-Mathijs Schoffelen and Robert Oostenveld
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

if isfield(volume, 'inside') && isfield(volume, 'outside')
  % reformat the inside field to a full volume so that it can be flipped/permuted as well
  tmp = [];
  tmp(volume.inside)  = 1;
  tmp(volume.outside) = 0;
  volume.inside  = tmp;
  volume.outside = [];
end

% determine the known volume parameters
param = parameterselection('all', volume);

% ensure that each of the volumes is 3D
for i=1:length(param)
  volume = setsubfield(volume, param{i}, reshape(getsubfield(volume, param{i}), volume.dim));
end

% permute the 3D volume so that the indices and head-coordinate axes correspond approximately
% the determination of the permutation order requires that voxel size differneces are accounted for
voxsiz = sqrt(sum(volume.transform(:,1:3).^2,1));

[dum, dim1] = max(abs(volume.transform(1,1:3))./voxsiz);
[dum, dim2] = max(abs(volume.transform(2,1:3))./voxsiz);
[dum, dim3] = max(abs(volume.transform(3,1:3))./voxsiz);
permutevec  = [dim1 dim2 dim3];
permutemat(4,4) = 1;
permutemat(permutevec(1),1) = 1;
permutemat(permutevec(2),2) = 1;
permutemat(permutevec(3),3) = 1;
transform = permutemat;
if length(unique(permutevec))<3
  ft_error('could not determine the correspondence between volume and headcoordinate axes');
else
  for i=1:length(param)
    volume = setsubfield(volume, param{i}, permute(getsubfield(volume, param{i}), permutevec));
  end
  volume.transform(:,1:3) = volume.transform(:,permutevec);
end

% update the dimensions of the volume
volume.dim = size(getsubfield(volume, param{1}));
if isfield(volume, 'xgrid')
  volume.xgrid = 1:volume.dim(1);
  volume.ygrid = 1:volume.dim(2);
  volume.zgrid = 1:volume.dim(3);
end

% subsequently flip the volume along each direction, so that the diagonal of the transformation matrix is positive
flipflags = zeros(1,3);
flipx = eye(4); flipx(1,1) = -1; flipx(1,4) = volume.dim(1)+1;
flipy = eye(4); flipy(2,2) = -1; flipy(2,4) = volume.dim(2)+1;
flipz = eye(4); flipz(3,3) = -1; flipz(3,4) = volume.dim(3)+1;
if volume.transform(1,1)<0
  flipflags(1) = 1;
  for i=1:length(param)
    volume = setsubfield(volume, param{i}, flip(getsubfield(volume, param{i}), 1));
  end
  volume.transform = volume.transform * flipx;
  transform        = transform * flipx;
end
if volume.transform(2,2)<0
  flipflags(2) = 1;
  for i=1:length(param)
    volume = setsubfield(volume, param{i}, flip(getsubfield(volume, param{i}), 2));
  end
  volume.transform = volume.transform * flipy;
  transform        = transform * flipy;
end
if volume.transform(3,3)<0
  flipflags(3) = 1;
  for i=1:length(param)
    volume = setsubfield(volume, param{i}, flip(getsubfield(volume, param{i}), 3));
  end
  volume.transform = volume.transform * flipz;
  transform        = transform * flipz;
end

if isfield(volume, 'inside') && isfield(volume, 'outside')
  % restore the inside and outside fields in their proper formats
  tmp = volume.inside(:);
  tmp(isnan(tmp)) = 0; % nan values also indicate outside the brain
  volume.inside   = find( tmp);
  volume.outside  = find(~tmp);
end

