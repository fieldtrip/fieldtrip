function err = hdf5write_safe(fname, name, val, options)

err = -1;
if  isempty(val)
    return;
end
if ~exist('options','var')
    options = '';
end

if iscell(val)
    HDF5_DatasetWriteStrings(fname, name, val)
else
    val = HDF5_Transpose(val, options);
    try
        if ~isempty(findstr('rw', options))
            HDF5_DatasetWrite(fname, name, val);
        else
            hdf5write(fname, name, val, 'WriteMode','append');
        end
    catch
        return;
    end
end
err = 0;