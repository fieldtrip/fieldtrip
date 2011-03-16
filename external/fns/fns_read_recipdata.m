function [data,compress,gridlocs,node_sizes,voxel_sizes] = fns_read_recipdata(infile)
% 
% This function reads the reciprocity data and the related
% information into the MATLAB workspace.

% $Copyright (C) 2010 by Hung Dang$


% Read in the reciprocity data
data = hdf5read(infile,'/recipdata/data');

% Read in other parameters
compress    = hdf5read(infile,'/recipdata/roicompress') + 1;
gridlocs    = hdf5read(infile,'/recipdata/gridlocs');
node_sizes  = hdf5read(infile,'/model/node_sizes');
voxel_sizes = hdf5read(infile,'/model/voxel_sizes');
