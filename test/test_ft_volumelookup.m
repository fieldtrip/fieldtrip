function test_ft_volumelookup

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_volumelookup ft_read_atlas atlas_lookup

atlasfilename = dccnpath('/home/common/matlab/fieldtrip/template/atlas/afni/TTatlas+tlrc.BRIK');
mrifilename   = dccnpath('/home/common/matlab/fieldtrip/external/spm8/templates/T1.nii');

mri = ft_read_mri(mrifilename);
mri.coordsys = 'mni';

cfg       = [];
cfg.atlas = atlasfilename;
cfg.roi   = {'Precentral Gyrus'};
mask1     = ft_volumelookup(cfg,mri);
cfg.roi   = {'Brodmann area 47'};
mask2     = ft_volumelookup(cfg,mri);
cfg.roi   = {'Precentral Gyrus';'Brodmann area 47'};
mask3     = ft_volumelookup(cfg,mri);

assert(isequal(sum(mask1(:)),8144));%8213));
assert(isequal(sum(mask2(:)),2096));%2093));
assert(isequal(sum(mask3(:)&(mask2(:)|mask1(:))),sum(mask3(:))));


% the following is just to check whether the functionality does not crash
% no idea who would use it anyway
cfg       = [];
cfg.roi   = [15 30 40;-40 -20 80];
cfg.sphere = [10 8];
mask4     = ft_volumelookup(cfg,mri);

mri.mask  = mask4;

cfg       = [];
cfg.atlas = atlasfilename;
cfg.maskparameter = 'mask';
mask5     = ft_volumelookup(cfg, mri);


atlasfilename = dccnpath('/home/common/matlab/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii');     
cfg       = [];
cfg.atlas = atlasfilename;
cfg.roi = 'Calcarine_R'; % right V1
mask6 = ft_volumelookup(cfg, mri);
assert(isequal(sum(mask6(:)),1861));

atlas_MNI = ft_read_atlas(atlasfilename);
atlas_MNI.coordsys = 'mni';
cfg = [];
cfg.roi = [52 -9 -45; 73 -37 -8];
cfg.atlas = atlasfilename;
cfg.output = 'single';
cfg.maxqueryrange = 29;
label_MNI = ft_volumelookup(cfg, atlas_MNI);
assert(size(label_MNI, 1) == 1);
cfg.output = 'label';
label_MNI = ft_volumelookup(cfg, atlas_MNI);
assert(size(label_MNI, 1) == 1);
cfg.output = 'multiple';
label_MNI = ft_volumelookup(cfg, atlas_MNI);
assert(size(label_MNI, 1) == 2);
