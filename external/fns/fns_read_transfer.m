function [trMatrix,elecStatus] = fns_read_transfer(datafile)
% Note that Data only stores the potential values for nodes, which  belong to
% the head

trMatrix      = hdf5read(datafile, '/forward/transfer_matrix');
elecStatus    = hdf5read(datafile, '/forward/status');
