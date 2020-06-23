function test_ft_read_and_plot_atlas

% WALLTIME 00:15:00
% MEM 3gb
% DEPENDENCY ft_read_atlas ft_sourceplot

% spm8 might have mexfile issues
ft_hastoolbox('spm12',1);
[ftver, ftpath] = ft_version;

% load MNI pial mesh
pial = load(fullfile(ftpath, '/template/anatomy/surface_pial_both.mat'));

% read and plot AFNI
atlas      = ft_read_atlas(fullfile(ftpath, '/template/atlas/afni/TTatlas+tlrc.HEAD'));
atlas.tissue = double(atlas.brick0(:,:,:,1));
atlas.dim = size(atlas.tissue);
atlas.coordsys = 'mni';

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, atlas);

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
cfg.method = 'surface';
ft_sourceplot(cfg, atlas, pial.mesh)

% read and plot AAL
atlas      = ft_read_atlas(fullfile(ftpath, '/template/atlas/aal/ROI_MNI_V4.nii'));

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, atlas);

cfg.method = 'surface';
ft_sourceplot(cfg, atlas, pial.mesh)

% read and plot Brainweb
brainweb = load(fullfile(ftpath, '/template/atlas/brainweb/brainweb_discrete.mat'));
atlas      = brainweb.atlas; clear brainweb

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, atlas);

cfg.method = 'surface';
ft_sourceplot(cfg, atlas, pial.mesh)

% read and plot JuBrain
atlas = ft_read_atlas(fullfile(ftpath, '/template/atlas/spm_anatomy/AllAreas_v18_MPM.mat'));

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, atlas);

cfg.method = 'surface';
ft_sourceplot(cfg, atlas, pial.mesh)

% read and plot VTPM
load(fullfile(ftpath, '/template/atlas/vtpm/vtpm.mat'));
atlas = vtpm; clear vtpm

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, atlas);

cfg.method = 'surface';
ft_sourceplot(cfg, atlas, pial.mesh)

% read and plot Brainnetome
atlas = ft_read_atlas(fullfile(ftpath, '/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii'));

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, atlas);

cfg.method = 'surface';
ft_sourceplot(cfg, atlas, pial.mesh)
