function rgn = hh_region_read(infile)
% hh_region_read - This routine reads all of the region information
% into a structure.
%
% Usage: rgn = hh_region_read(infile)
% where 
%       1. infile:    An input HDF5 data file
%       2. rgn:       A returned region 
% 
% $Author: Hung Dang$
% $Id: hh_region_read.m$
% $Log$
% Revision 1.1 Mon Aug  9 20:27:29 MDT 2010, hungptit
% First update

% $$$ Init the region data structure
rgn = struct('locations',[],'gridlocs',[],'values',[],'status',[], ...
             'node_sizes',[],'voxel_sizes',[]);

% $$$ Read in the region information
rgn.locations   = hdf5read(infile,'/region/locations');
rgn.gridlocs    = hdf5read(infile,'/region/gridlocs');
rgn.values      = hdf5read(infile,'/region/values');
rgn.status      = hdf5read(infile,'/region/status');
rgn.node_sizes  = hdf5read(infile,'/region/node_sizes');
rgn.voxel_sizes = hdf5read(infile,'/region/voxel_sizes');
rgn.info        = ''; 
% $$$ rgn.info        = hdf5read(infile,'/region/info');