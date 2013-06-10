function test_bug2193

% TEST test_bug2193
% TEST ft_read_atlas

filename = dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2193/aal_MNI_V4.nii');

aal = ft_read_atlas(filename);

assert(all(~cellfun(@isempty, aal.tissuelabel)), 'there is an empty tissuelabel');
assert(max(aal.tissue(:))==length(aal.tissuelabel), 'inconsistent number of tissues');



