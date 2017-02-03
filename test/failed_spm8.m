function failed_spm8

% WALLTIME 00:10:00
% MEM 3gb

% TEST test_spm8 ft_volumenormalise ft_volumesegment ft_volumedownsample mni2tal tal2mni

mrifile = dccnpath('/home/common/matlab/fieldtrip/data/bauer_m.mri');
mri     = ft_read_mri(mrifile);

[ftver, ftpath] = ft_version;

% start with a clean path
rmpath(fullfile(ftpath, 'external', 'spm2'));
rmpath(fullfile(ftpath, 'external', 'spm8'));
rmpath(fullfile(ftpath, 'external', 'spm12'));

%% ft_volumenormalise

cfg            = [];
cfg.nonlinear  = 'no';

clear fun;
cfg.spmversion = 'spm2';
n2 = ft_volumenormalise(cfg, mri);
rmpath(spm('dir'));

clear fun;
cfg.spmversion = 'spm8';
n8 = ft_volumenormalise(cfg, mri);
rmpath(spm('dir'));


%% ft_volumesegment

clear fun;
cfg             = [];
cfg.spmversion  = 'spm2';
s2 = ft_volumesegment(cfg, mri);
rmpath(spm('dir'));

clear fun;
cfg.spmversion  = 'spm8';
s8 = ft_volumesegment(cfg, mri);
rmpath(spm('dir'));


%% ft_volumerealign

clear fun;
cfg               = [];
cfg.method        = 'spm';
cfg.spmversion    = 'spm2';
r2 = ft_volumerealign(cfg, mri, mri);
rmpath(spm('dir'));

clear fun;
cfg.spmversion    = 'spm8';
r8 = ft_volumerealign(cfg, mri, mri);
rmpath(spm('dir'));


%% ft_volumedownsample
tmp            = randn(256,256,256);
mri.pow        = tmp;
mri.pow(1)     = 0;

cfg            = [];
cfg.downsample = 2;

clear fun;
cfg.spmversion = 'spm2';
cfg.smooth     = 'no';
d2             = ft_volumedownsample(cfg, mri);
cfg.smooth     = 5;
d2s            = ft_volumedownsample(cfg, mri);
rmpath(spm('dir'));

clear fun;
cfg.spmversion = 'spm8';
cfg.smooth     = 'no';
d8             = ft_volumedownsample(cfg, mri);
cfg.smooth     = 5;
d8s            = ft_volumedownsample(cfg, mri);
rmpath(spm('dir'));


%% mni2tal and tal2mni

[ftver, ftpath] = ft_version;
cd(fullfile(ftpath, 'private'));

inpoints  = randn(100,3);
outpoints = mni2tal(inpoints);
rmpath(spm('dir'));

inpoint2s = tal2mni(outpoints);
rmpath(spm('dir'));

%% prepare_dipole_grid

%% prepare_mesh_segmentation
