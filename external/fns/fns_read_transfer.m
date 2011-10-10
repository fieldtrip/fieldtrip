function forwardSolutions = fns_read_transfer(datafile)
forwardSolutions = hdf5read(datafile, '/forward/transfer_matrix');
% forwardSolutions = struct('Data', [], 'Compress', [], ...
%     'ModelSizes', [], 'VoxelSizes', []);
% 
% % Note that Data only stores the potential values for nodes, which are belong to
% % the head.
% forwardSolutions.Data       = hdf5read(datafile, '/forward/transfer_matrix');
% 
% % These parameters are used to compute the lead field matrices.
% forwardSolutions.Compress   = hdf5read(datafile, '/sparse/compress');
% forwardSolutions.ModelSizes = hdf5read(datafile, '/model/node_sizes');
% forwardSolutions.VoxelSizes = hdf5read(datafile, '/model/voxel_sizes');
end
