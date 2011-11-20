function [trMatrix,elecStatus] = fns_read_transfer(datafile)
% Note that Data only stores the potential values for nodes, which  belong to
% the head

trMatrix      = hdf5read(datafile, '/forward/transfer_matrix');
elecStatus    = hdf5read(datafile, '/forward/status');

% % These parameters are used to compute the lead field matrices.
% forwardSolutions.Compress   = hdf5read(datafile, '/sparse/compress');
% forwardSolutions.ModelSizes = hdf5read(datafile, '/model/node_sizes');
% forwardSolutions.VoxelSizes = hdf5read(datafile, '/model/voxel_sizes');
end
