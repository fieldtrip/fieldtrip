function test_bug2193

% TEST test_bug2193
% TEST ft_read_atlas

%%%%%% from http://fmri.wfubmc.edu

filename = dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2193/wfu/aal_MNI_V4.nii');

aal = ft_read_atlas(filename);
assert(all(~cellfun(@isempty, aal.tissuelabel)), 'there is an empty tissuelabel');
assert(max(aal.tissue(:))==length(aal.tissuelabel), 'inconsistent number of tissues');

cfg = [];
cfg.atlas = filename;
bbl = ft_prepare_atlas(cfg);


%%%%%% from http://www.cyceron.fr
filename = dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2193/aal/ROI_MNI_V4.nii');

aal = ft_read_atlas(filename);
assert(all(~cellfun(@isempty, aal.tissuelabel)), 'there is an empty tissuelabel');
assert(max(aal.tissue(:))==length(aal.tissuelabel), 'inconsistent number of tissues');

cfg = [];
cfg.atlas = filename;
bbl = ft_prepare_atlas(cfg);
