function err = HDF5_DatasetWrite(fname, location, data)

if exist(fname,'file')
    fid = H5F.open(fname, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
else
    fid = H5F.create(fname, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
end
if fid<0
    err = -1;
    return;
end

% Create group where dataset is located
group = fileparts(location);
group(group=='\') = '/';
gid = CreateGroup(fid, group);
if gid < 0
    err = -1;
    return;
end

% Update dataset
err = WriteDataset(fid, location, data);




% ----------------------------------------------------------------
function gid = CreateGroup(fid, group)
if group == '/'
    gid = H5G.open(fid, '/');
    return;
end

try
    gid = H5G.open(fid, group);
catch
    [rootgroup, group] = fileparts(group);
    rootgroup(rootgroup=='\') = '/';
    gid = CreateGroup(fid, rootgroup);
    gid = H5G.create(gid, group, 'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
end



% --------------------------------------------------------------
function err = WriteDataset(fid, location, data)

err = 0;

% data = HDF5_Transpose(data);
dims = size(data);
h5_dims = fliplr(dims);
unlimited = H5ML.get_constant_value('H5S_UNLIMITED');
h5_maxdims = [unlimited, unlimited];

if isnumeric(data)
    type_id = H5T.copy('H5T_NATIVE_DOUBLE');
elseif ischar(data)
    type_id = H5T.copy('H5T_C_S1'); 
else
    err = -1;
    return;
end

space_id = H5S.create_simple(2, h5_dims, h5_maxdims);
proplist_id = H5P.create('H5P_DATASET_CREATE');

% Open or create dataset
try
    dset_id = H5D.open(fid, location);
    H5D.set_extent(dset_id, h5_dims);
catch
    H5P.set_chunk(proplist_id, h5_dims);
    dset_id = H5D.create(fid, location, type_id, space_id, proplist_id);
end
file_space_id = H5D.get_space(dset_id);

% Write dataset
% H5D.write(dset_id, type_id, space_id, file_space_id, proplist_id, data);
H5D.write(dset_id, type_id, space_id, file_space_id, 'H5P_DEFAULT', data);

% Close space, type, dataset and file
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);
H5F.close(fid);



% ----------------------------------------------------------------
function [dset_id, type_id, space_id, proplist_id] = CreateDataset(fid, location, h5_dims)

h5_maxdims = h5_dims;

type_id = H5T.copy('H5T_NATIVE_DOUBLE');
space_id = H5S.create_simple(2, h5_dims, h5_maxdims);
proplist_id = H5P.create('H5P_DATASET_CREATE');
H5P.set_chunk(proplist_id, [h5_dims(1), h5_dims(2)]);

dset_id = H5D.create(fid, location, type_id, space_id, proplist_id);


