function [transform] = pos2transform(pos)

% POS2TRANSFORM reconstructs a transformation matrix from an ordered list 
% of positions.
%
% Use as
%  [transform] = pos2transform(pos, dim) where pos is an ordered list of positions
%  and should specify a full 3D volume 
%
% The output transform is a 4x4 homogenous transformation matrix which transforms
% from 'voxelspace' into the positions provided in the input

% Copyright (C) 2009, Jan-Mathijs Schoffelen

dim = pos2dim3d(pos);
x   = 1:dim(1);
y   = 1:dim(2);
z   = 1:dim(3);
[X,Y,Z] = ndgrid(x, y, z);
ind = [X(:) Y(:) Z(:)];
ind = ind';ind(4,:) = 1;
pos = pos';pos(4,:) = 1;
transform = pos/ind;
