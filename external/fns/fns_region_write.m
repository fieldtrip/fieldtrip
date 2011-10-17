%% TODO: Need to remove this routine.
function fns_region_write(rgnfile,rgn)
% 
% This routine write the region information to the HDF5 file
%
% rgnfile The output region file.
% rgn The region of interest data structure.

% @author Hung Dang, July 13, 2010

% Write the description of stored data to the output file
hdf5write(rgnfile,'/region/info',rgn.info);

% Write the locations matrix to the /region/locations dataset
hdf5write(rgnfile,'/region/locations',rgn.locations,'WriteMode', ...
          'append');

% Write the gridlocs matrix  to the /region/gridlocs dataset
hdf5write(rgnfile,'/region/gridlocs',rgn.gridlocs,'WriteMode', ...
          'append');

% Write the voxel_sizes vector to the /region/voxel_sizes dataset
hdf5write(rgnfile,'/region/voxel_sizes',rgn.voxel_sizes, ...
          'WriteMode','append');

% Write the node_sizes vector  to the /region/node_sizes dataset
hdf5write(rgnfile,'/region/node_sizes',rgn.node_sizes,'WriteMode', ...
          'append');

% Write the status vector  to the /region/status dataset
hdf5write(rgnfile,'/region/status',rgn.status,'WriteMode', ...
          'append');

% Write the values vector  to the /region/status dataset
hdf5write(rgnfile,'/region/values',rgn.values,'WriteMode', ...
          'append');
