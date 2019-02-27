function [transform] = pos2transform(pos, dim)

% POS2TRANSFORM reconstructs a transformation matrix from an ordered list 
% of positions.
%
% Use as
%   [transform] = pos2transform(pos, dim)
% where pos is an ordered list of positions that should specify a full 3D volume.
%
% The output transform is a 4x4 homogenous transformation matrix which transforms
% from 'voxelspace' into the positions provided in the input
%
% See also POS2DIM

% Copyright (C) 2009, Jan-Mathijs Schoffelen

if nargin>1
  % do nothing
else
  dim = pos2dim(pos);
end
x   = 1:dim(1);
y   = 1:dim(2);
z   = 1:dim(3);
[X,Y,Z] = ndgrid(x, y, z);
ind = [X(:) Y(:) Z(:)];
ind = ind'; ind(4,:) = 1;
pos = pos'; pos(4,:) = 1;

% build in some robustness against nans
sel = sum(isfinite(pos))==4;

transform = pos(:,sel)/ind(:,sel);

