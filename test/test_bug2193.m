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


% The following section depends on the current (10 June 2013) situation with spm8
% on home/common, which will probably not persist forever. Therefore this section
% should not be included in the automatic regression test.

if false
  cfg.atlas = '/home/common/matlab/spm8/toolbox/wfu_pickatlas_3.03-old2/MNI_atlas_templates/aal_MNI_V4.nii'
  aal1a = ft_prepare_atlas(cfg);
  aal1b = ft_read_atlas(cfg.atlas);
  
  cfg.atlas = '/home/common/matlab/spm8/toolbox/wfu_pickatlas/MNI_atlas_templates/aal_MNI_V4.nii'
  aal2a = ft_prepare_atlas(cfg);
  aal2b = ft_read_atlas(cfg.atlas);
  
  cfg.atlas = '/home/common/matlab/spm8/toolbox/wfu_pickatlas-old/MNI_atlas_templates/aal_MNI_V4.img'
  aal3a = ft_prepare_atlas(cfg);
  aal3b = ft_read_atlas(cfg.atlas);
end