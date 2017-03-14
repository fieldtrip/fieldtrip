function test_ft_volumelookup

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_volumelookup
% TEST ft_read_atlas
% TEST atlas_lookup

atlasfilename = dccnpath('/home/common/matlab/fieldtrip/template/atlas/afni/TTatlas+tlrc.BRIK');
mrifilename   = dccnpath('/home/common/matlab/fieldtrip/external/spm8/templates/T1.nii');

mri = ft_read_mri(mrifilename);

cfg       = [];
cfg.atlas = atlasfilename;
cfg.inputcoord = 'mni';
cfg.roi   = {'Precentral Gyrus'};
mask1     = ft_volumelookup(cfg,mri);
cfg.roi   = {'Brodmann area 47'};
mask2     = ft_volumelookup(cfg,mri);
cfg.roi   = {'Precentral Gyrus';'Brodmann area 47'};
mask3     = ft_volumelookup(cfg,mri);

assert(isequal(sum(mask1(:)),8214));
assert(isequal(sum(mask2(:)),2102));
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
cfg.inputcoord = 'mni';
cfg.maskparameter = 'mask';
mask5     = ft_volumelookup(cfg, mri);


atlasfilename = dccnpath('/home/common/matlab/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii');     
cfg       = [];
cfg.atlas = atlasfilename;
cfg.inputcoord = 'mni';
cfg.roi = 'Calcarine_R'; % right V1
mask6 = ft_volumelookup(cfg, mri);
assert(isequal(sum(mask6(:)),1861));
