function [voxel, head] = cornerpoints(dim, transform)

% CORNERPOINTS returns the eight corner points of an anatomical volume
% in voxel and in head coordinates
%
% Use as
%   [voxel, head] = cornerpoints(dim, transform)
% which will return two 8x3 matrices.

% determine the corner points of the volume in voxel space
voxel = [
  1      1      1
  dim(1)     1      1
  dim(1) dim(2)     1
  1  dim(2)     1
  1      1  dim(3)
  dim(1)     1  dim(3)
  dim(1) dim(2) dim(3)
  1  dim(2) dim(3)
  ];

% determine the corner points of the volume in plotting space
head = warp_apply(transform, voxel);
