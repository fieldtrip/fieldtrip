function test_bug1573

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_read_mri ft_volumnrealign ft_sourceinterpolate

template = dccnpath('/home/common/matlab/spm8/canonical/single_subj_T1.nii');
template_mri = ft_read_mri(template);

cfg = [];
cfg.method = 'fiducial';
cfg.fiducial.nas = [45 106 17];
cfg.fiducial.lpa = [88 49 11];
cfg.fiducial.rpa = [3 49 11];
mri = ft_volumerealign(cfg,template_mri);

%%
return
% these files johzum has but are not on svn.
load('gava_sourcenorm/grandAVG_SourceDiff_OHNEdysphagie__beta.mat');
% below requires too much memory

cfg = [];
cfg.downsample = 10;
cfg.parameter = 'avg.avg.pow';
sourceDiff_int = ft_sourceinterpolate(cfg,grandavg,mri);

cfg = [];
cfg.method = 'ortho';
cfg.interactive ='yes';
cfg.funparameter ='avg.avg.pow';
figure; ft_sourceplot(cfg,sourceDiff_int);


