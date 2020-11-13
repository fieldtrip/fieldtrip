function [gid, fid] = HDF5_GroupOpen(fileobj, location)

% Find data location
if ischar(fileobj)
    fid = H5F.open(fileobj);
else
    fid = fileobj;
end
try
    gid = H5G.open(fid, location);
catch
    gid.double = -1;
end


