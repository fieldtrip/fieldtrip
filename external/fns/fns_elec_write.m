function fns_elec_write(pnt, vsize, dimpos, elecfile)
%
% This function write the electrodes locations on disk and other
% information
%   1. pnt is the location of the electrodes [NX3]
%   3. voxel_sizes is the dimension of a voxel in the cartesian coordinates in mm [3X1]
%   4. nodes_sizes is the dimension of the dipoles grid [NX NY NZ]

% $Copyright (C) 2010 by Hung Dang$

hdf5write(elecfile, '/electrodes/locations', pnt);

hdf5write(elecfile, '/electrodes/gridlocs', int32(pnt), ...
          'WriteMode', 'append');       % Assume the electrodes
                                        % locations are mm and the
                                        % voxel sizes is 1mm x 1mm
                                        % x 1mm.

hdf5write(elecfile, '/electrodes/voxel_sizes', vsize, ...
          'WriteMode','append');

hdf5write(elecfile, '/electrodes/node_sizes', int32(dimpos), ...
          'WriteMode', 'append');

