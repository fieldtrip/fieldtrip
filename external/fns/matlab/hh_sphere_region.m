function newrgn = hh_sphere_region(rgn,voxel_sizes,r,center)
% hh_sphere_region - This function return the subset of the given
% region which lie inside the sphere center at c with radius r
% 
% Usage: newrgn = hh_sphere_region(rgn,r,c)
% 
% $ID: hh_sphere_region.m$
% $AUTHOR: Copyright (C) by Hung Dang$
% $DATE: Thu Aug  5 14:29:18 MDT 2010$
%
% %LOG%
% Revision 1.1 Mon Sep 27 17:54:50 MDT 2010, hungptit
% Update the documentation 
% 

% Update parameters
X = double(rgn);

% Compute the distance from the center to all nodes
dx = (X(:,1) - center(1)) * voxel_sizes(1);
dy = (X(:,2) - center(2)) * voxel_sizes(2);
dz = (X(:,3) - center(3)) * voxel_sizes(3);
dis = dx.^2 + dy.^2 + dz.^2;
idx = dis < r^2;

% Update the new region
newrgn = rgn(idx,:);
