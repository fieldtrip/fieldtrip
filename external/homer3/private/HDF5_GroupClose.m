function HDF5_GroupClose(fileobj, gid, fid)

H5G.close(gid);
if ischar(fileobj)
    H5F.close(fid);
end
