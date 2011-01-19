% Syntax - fwhm = hh_fwhm(cost,decompress,node_sizes,voxel_sizes) 
%
% $AUTHOR: Copyright (C) by Hung Dang$
% $DATE: Tue Sep 14 18:56:22 MDT 2010$
% $ID: hh_fwhm.m$
% 
% 
% $LOG$
% Revision 1.1 Tue Sep 14 18:56:54 MDT 2010, hungptit
% First update and this routine has been checked with the C code.

function fwhm = hh_fwhm(cost,decompress,node_sizes,voxel_sizes)

% Update parameters
ROW = node_sizes(1);
COL = node_sizes(2);
dx = voxel_sizes(1);
dy = voxel_sizes(2);
dz = voxel_sizes(3);
threshold = max(cost) * 0.5;

% Find max node only
idx = decompress(cost == max(cost));
row_max = mod(idx,ROW);
idx = (idx - row_max) / ROW;
col_max = mod(idx,COL);
slice_max = (idx - col_max) / COL;

% Find all valid nodes
idx = decompress(cost > threshold);
row = mod(idx,ROW);
idx = (idx - row) / ROW;
col = mod(idx,COL);
slice = (idx - col) / COL;

% Compute the distance
lx = double(row_max - row) * dx;
ly = double(col_max - col) * dy;
lz = double(slice_max - slice) * dz;
fwhm = sqrt(max(lx .* lx + ly .* ly + lz .* lz));