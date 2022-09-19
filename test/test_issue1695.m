function test_issue1695

% MEM 2gb
% WALLTIME 01:00:00
% DEPENDENCY ft_version ft_read_mri ft_volumereslice ft_volumenormalise ft_warp_apply spm_file spm_preproc_run

%% Init
[ftver, ftpath] = ft_version;

%% Read sourcemodel
dat = load(fullfile(ftpath, 'template/sourcemodel/standard_sourcemodel3d10mm.mat'));
sourcemodel = dat.sourcemodel;

%% Normalise #1 - FieldTrip
mri = ft_read_mri(fullfile(ftpath, 'template/anatomy/single_subj_T1.nii'),'dataformat','nifti_spm');
mri.coordsys = 'ras';
cfg     = [];
cfg.dim = mri.dim;
mri     = ft_volumereslice(cfg,mri);

% this seemed to initially crash
cfg = [];
cfg.spmversion = 'spm12';
cfg.spmmethod = 'new';
cfg.nonlinear = 'no';
normalise = ft_volumenormalise(cfg, mri);
pos       = ft_warp_apply(normalise.params, sourcemodel.pos, 'sn2individual');

% this worked before
cfg.nonlinear = 'yes';
normalise2 = ft_volumenormalise(cfg, mri);
pos2       = ft_warp_apply(normalise2.params, sourcemodel.pos, 'sn2individual');
