function [ data ] = read_ah5_data( filename, hdr, begsample, endsample, chanindx)

offset = (begsample-1)*hdr.nChans;
nsamples = (endsample-begsample+1);
block = nsamples * hdr.nChans;

fid = H5F.open(filename);
dset_id = H5D.open(fid, '/Blocks/1/data');
mem_space = H5S.create_simple(1, block, block);
space = H5D.get_space(dset_id);

H5S.select_hyperslab(space, 'H5S_SELECT_SET', offset, [], [], block);
data = H5D.read(dset_id, 'H5ML_DEFAULT', mem_space, space, 'H5P_DEFAULT');
data = reshape(data, hdr.nChans, nsamples);
H5D.close(dset_id);
H5F.close(fid);

end

