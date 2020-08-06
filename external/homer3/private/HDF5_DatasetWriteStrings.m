function HDF5_DatasetWriteStrings(fname, location, data)

if ~exist(fname, 'file')
    fid = H5F.create(fname, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
else
    fid = H5F.open(fname,'H5F_ACC_RDWR','H5P_DEFAULT');
end
data = cell2str_new(data);

filetype = H5T.copy('H5T_FORTRAN_S1'); 
H5T.set_size(filetype, size(data,2)); 
memtype = H5T.copy('H5T_C_S1'); 
H5T.set_size(memtype, size(data,2)); 

% Create dataspace. Setting maximum size to [] sets the maximum 
% size to be the current size. 
space = H5S.create_simple(1, size(data,1), []); 

% Create the dataset and write the string data to it. 
dset = H5D.create(fid, location, filetype, space, 'H5P_DEFAULT'); 

% Transpose the data to match the layout in the H5 file to match C 
% generated H5 file. 
H5D.write(dset, memtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT',  data'); 

% Close and release resources. 
H5D.close(dset); 
H5S.close(space); 
H5T.close(filetype); 
H5T.close(memtype); 
H5F.close(fid); 


