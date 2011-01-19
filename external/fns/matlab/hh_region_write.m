function hh_region_write(rgnfile,rgn)
% @brief This routine write the region information to the HDF5 file
% @author Hung Dang
% @date July 13, 2010
%
% @param rgnfile The output region file.
% @param rgn The region of interest data structure.
% @param description The description about the stored data.
% $Log$
% Revision 1.1 Mon Aug  9 20:22:21 MDT 2010, hungptit
% Update the code to write the values to the region file
% Revision 1.2 Mon Aug  9 20:31:22 MDT 2010, hungptit
% Update code to write the region information according to the new
% region data structure.

% $$$ Write the description of stored data to the output file
hdf5write(rgnfile,'/region/info',rgn.info);

% $$$ Write the locations matrix to the /region/locations dataset
hdf5write(rgnfile,'/region/locations',rgn.locations,'WriteMode', ...
          'append');

% $$$ Write the gridlocs matrix  to the /region/gridlocs dataset
hdf5write(rgnfile,'/region/gridlocs',rgn.gridlocs,'WriteMode', ...
          'append');

% $$$ Write the voxel_sizes vector to the /region/voxel_sizes dataset
hdf5write(rgnfile,'/region/voxel_sizes',rgn.voxel_sizes, ...
          'WriteMode','append');

% $$$ Write the node_sizes vector  to the /region/node_sizes dataset
hdf5write(rgnfile,'/region/node_sizes',rgn.node_sizes,'WriteMode', ...
          'append');

% $$$ Write the status vector  to the /region/status dataset
hdf5write(rgnfile,'/region/status',rgn.status,'WriteMode', ...
          'append');

% $$$ Write the values vector  to the /region/status dataset
hdf5write(rgnfile,'/region/values',rgn.values,'WriteMode', ...
          'append');
